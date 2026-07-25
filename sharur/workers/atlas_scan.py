"""Atlas scan worker executor.

Implements the eight-step production worker contract for `atlas_genome_read`
tasks:

1. register a distinct agent identity and exact ``profile:<name>`` capability
2. claim compatible tasks over Ops HTTP
3. launch the resolved provider/model/effort from the pinned review policy
4. pass only the frozen bounded, sequence-free task input
5. heartbeat and checkpoint under the attempt fence
6. validate and submit the required structured record
7. complete the task with the same attempt token
8. preserve stderr, exit status, model metadata, prompt hash, result hash

One genome per task. The worker walks that genome's adaptive packet frames to
`complete`, accumulates typed candidate occurrences, writes the exact coverage
manifest, then submits candidates before exactly one unit disposition — the
order the review-output contract requires.
"""

from __future__ import annotations

import hashlib
import json
import logging
import random
import signal
import time
from pathlib import Path
from typing import Any

import requests

from sharur.atlas import frame_coverage_receipt, write_genome_coverage_manifest
from sharur.ops.client import SharurOps
from sharur.query.client import SharurQuery
from sharur.review.models import load_review_policy
from sharur.workers.model_cli import ModelError, ModelRateLimited, ModelRun, run_profile

LOGGER = logging.getLogger("sharur.workers.atlas_scan")

TASK_TYPE = "atlas_genome_read"
CANDIDATE_TYPE_PREFIX = "atlas-scan"
SIGNATURE_SCHEMA = "atlas-scan-signature/1.0"

# What the model is allowed to emit. Enforced server-side as a schema for the
# OpenAI driver; stated in-prompt for the Anthropic driver.
SCAN_OUTPUT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["candidates"],
    "properties": {
        "candidates": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": [
                    "candidate_type",
                    "signature",
                    "evidence",
                    "subject_refs",
                ],
                "properties": {
                    "candidate_type": {"type": "string"},
                    "signature": {"type": "object"},
                    "evidence": {"type": "object"},
                    "subject_refs": {"type": "object"},
                    "reason_codes": {"type": "array", "items": {"type": "string"}},
                    "uncertainty": {"type": "object"},
                },
            },
        },
        "notes": {"type": "string"},
    },
}

SCAN_SYSTEM_PROMPT = """\
You are an Atlas scanner reading one bounded frame of one microbial genome.

You will receive a JSON packet containing contigs from exactly one genome, each
with its proteins. Every protein carries: protein_id, gene_index, start, end,
strand, length_aa, gc_content, observed_annotations, predicates, named_calls,
and loci. There are no biological sequences and you must never ask for them.

Your job is to report NOTABLE architecture as typed candidate occurrences. Do
not summarize the genome. Do not report routine housekeeping. Emit a candidate
only when the neighborhood shows something worth a reviewer's time: an unusual
gene cluster, a complete pathway whose parts are co-located, an unexpected
system in this lineage, a striking multi-domain architecture, or a locus of
co-directional uncharacterized proteins flanked by something informative.

SCIENTIFIC CONTRACT — these are hard rules:

1. OBSERVED evidence and NAMED calls are different things. An HMM hit in
   observed_annotations is an observation. A functional name requires either a
   purpose-built caller's output (named_calls) or co-annotation plus
   neighborhood support that you state explicitly.
2. Known annotation confounds. Check these before naming anything:
   - hyddb group 4 / energy-converting / H2-evolving NiFe calls are usually
     respiratory Complex I (NuoD/NqoD homology). Require NiFeSe_Hases on the
     same protein; Complex1_49kDa or Complex1_30kDa means Complex I, not a
     hydrogenase.
   - Any annotation flagged hyddb_needs_curation is low-confidence by the
     caller's own admission. Never emit a named hydrogenase class from one.
   - RuBisCO_large alone is not carbon fixation. Form IV RuBisCO-like proteins
     do methionine salvage. Require PRK in the same genome; RuBisCO_small
     further discriminates form I.
   - Ald_Xan_dh_C contains CoxL but also xanthine dehydrogenase and aldehyde
     oxidoreductase. Never call CO oxidation from the family alone.
   - Any accession appearing at >10 hits per genome describes a fold, not a
     function.
3. Phrase domain-compatible interpretations as hypotheses and mark them
   UNVERIFIED. Cite stable identifiers for everything.
4. Absence in a MAG is not absence in the organism. Do not report negatives.
5. Keep signatures typed, structured, and free of interpretive prose. Prose
   belongs in evidence.

OUTPUT

Return one JSON object with a "candidates" array. Each candidate needs:
  candidate_type   short kebab-case class, e.g. "co-located-pathway",
                   "novel-gene-cluster", "multidomain-architecture",
                   "mobile-element-cargo", "defense-locus"
  signature        typed structured fields ONLY (identifiers, accessions,
                   counts, coordinates). This is the exact reduction key —
                   two occurrences of the same finding must produce byte-
                   identical signatures. No prose, no free text.
  evidence         your reasoning, the observations supporting it, and any
                   UNVERIFIED hypotheses
  subject_refs     {"protein_ids": [...], "contig_ids": [...]}
  reason_codes     optional short tags
  uncertainty      optional {"confidence": "low|medium|high", "caveats": [...]}

An empty candidates array is a valid and common answer. Most frames are
housekeeping. Report nothing rather than inflating a routine frame.
"""


class LeaseLost(RuntimeError):
    """Another attempt owns this task; abandon it immediately."""


def _stable_key(*parts: str) -> str:
    return hashlib.sha256("\x1f".join(parts).encode("utf-8")).hexdigest()[:32]


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _assert_sequence_free(payload: dict[str, Any], genome_id: str) -> None:
    """Defense in depth. The operator and census both guarantee this already."""
    for contig in payload.get("contigs", []) or []:
        for protein in contig.get("proteins", []) or []:
            if isinstance(protein, dict) and "sequence" in protein:
                raise RuntimeError(
                    f"packet for {genome_id} carried sequence data; refusing to send to a model"
                )


def _is_conflict(exc: Exception) -> bool:
    resp = getattr(exc, "response", None)
    return resp is not None and resp.status_code == 409


class AtlasScanWorker:
    def __init__(
        self,
        *,
        ops_url: str,
        query_url: str,
        agent_id: str,
        profile: str,
        policy_path: str | None = None,
        campaign_id: str | None = None,
        ops_token: str | None = None,
        query_token: str | None = None,
        lease_seconds: int = 900,
        model_timeout: int = 1800,
        max_frames: int | None = None,
        dry_run: bool = False,
    ) -> None:
        self.agent_id = agent_id
        self.profile_name = profile
        self.campaign_id = campaign_id
        self.lease_seconds = lease_seconds
        self.model_timeout = model_timeout
        self.max_frames = max_frames
        self.dry_run = dry_run

        self.policy = load_review_policy(Path(policy_path) if policy_path else None)
        if profile not in self.policy.profiles:
            raise SystemExit(
                f"profile {profile!r} not in policy {self.policy.name!r}; "
                f"available: {sorted(self.policy.profiles)}"
            )
        self.profile = self.policy.profiles[profile]

        self.ops = SharurOps(ops_url, agent_id=agent_id, api_token=ops_token)
        self.query = SharurQuery(query_url, api_token=query_token)
        self._stop = False

    # ---------------------------------------------------------------- lifecycle

    def install_signal_handlers(self) -> None:
        def _handle(signum, _frame):
            LOGGER.info("signal %s received; finishing current task then stopping", signum)
            self._stop = True

        signal.signal(signal.SIGTERM, _handle)
        signal.signal(signal.SIGINT, _handle)

    def register(self) -> None:
        """Step 1 — identity plus the exact profile capability."""
        caps = ["atlas_reader", self.profile.capability]
        self.ops.register_agent(
            self.agent_id,
            role="worker",
            capabilities=caps,
            max_concurrent_tasks=1,
            capacity_cpu_slots=1,
        )
        LOGGER.info(
            "registered %s with capabilities %s -> %s/%s effort=%s",
            self.agent_id,
            caps,
            self.profile.provider,
            self.profile.model,
            self.profile.reasoning_effort,
        )

    def run_forever(self, *, idle_sleep: float = 5.0, max_tasks: int | None = None) -> int:
        self.register()
        done = 0
        rate_limit_backoff = 0.0
        while not self._stop:
            if max_tasks is not None and done >= max_tasks:
                break
            try:
                task = self.ops.claim_next_task(
                    campaign_id=self.campaign_id,
                    task_types=[TASK_TYPE],
                    lease_seconds=self.lease_seconds,
                )
            except requests.HTTPError as exc:
                LOGGER.warning("claim failed: %s", exc)
                time.sleep(idle_sleep)
                continue

            if task is None:
                time.sleep(idle_sleep)
                continue

            try:
                self.run_task(task)
                done += 1
                rate_limit_backoff = 0.0
            except LeaseLost:
                LOGGER.warning("lease lost on task %s; abandoning", task.get("id"))
            except ModelRateLimited as exc:
                # The subscription window is exhausted. This is the steady
                # state, not an error: hold nothing, release the task for
                # retry, and back off until the window rolls.
                rate_limit_backoff = min(max(rate_limit_backoff * 2, 60.0), 1800.0)
                jitter = random.uniform(0, rate_limit_backoff * 0.1)
                LOGGER.warning(
                    "rate limited (%s); releasing task %s and sleeping %.0fs",
                    exc,
                    task.get("id"),
                    rate_limit_backoff + jitter,
                )
                self._release(task, str(exc), retry_delay=int(rate_limit_backoff))
                time.sleep(rate_limit_backoff + jitter)
            except Exception as exc:  # noqa: BLE001 - worker must not die on one genome
                LOGGER.exception("task %s failed", task.get("id"))
                self._release(task, f"{type(exc).__name__}: {exc}", retry_delay=30)
        return done

    def _release(self, task: dict[str, Any], error: str, *, retry_delay: int) -> None:
        try:
            self.ops.fail_task(
                task["id"],
                error=error[:2000],
                retryable=True,
                retry_delay_seconds=retry_delay,
            )
        except Exception:  # noqa: BLE001
            LOGGER.warning("could not release task %s", task.get("id"), exc_info=True)

    # -------------------------------------------------------------------- task

    def run_task(self, task: dict[str, Any]) -> None:
        params = task.get("params") or {}
        task_id = task["id"]
        genome_id = params["genome_id"]
        unit_id = params["unit_id"]
        dataset_id = params["dataset_id"]
        campaign_id = task.get("campaign_id") or self.campaign_id
        checkpoint_key = params.get("checkpoint_key", "atlas_progress")
        checkpoint_interval = int(params.get("checkpoint_interval_frames", 1) or 1)

        LOGGER.info("task %s genome %s (%s proteins)", task_id, genome_id, params.get("n_proteins"))

        cursor, frames, candidates, contig_state = self._resume(task_id, checkpoint_key)

        frame_index = len(frames)
        while True:
            if self.max_frames is not None and frame_index >= self.max_frames:
                LOGGER.info("max_frames=%s reached; stopping early", self.max_frames)
                break

            envelope = self.query.genome_packet(
                genome_id,
                cursor=cursor,
                max_model_payload_bytes=int(params["packet_model_payload_bytes"]),
                all_annotations=bool(params.get("packet_all_annotations", True)),
            )
            raw = envelope["raw"]

            # Step 4 — the input must be exactly the frozen contract.
            if raw.get("dataset_id") != dataset_id:
                raise RuntimeError(f"packet dataset_id drift on {genome_id}")
            if raw.get("bin_id") != genome_id:
                raise RuntimeError(f"packet bin_id drift on {genome_id}")
            expected_hash = params.get("packet_packing_contract_hash")
            if expected_hash and raw.get("packing_contract_hash") != expected_hash:
                raise RuntimeError(f"packing contract hash drift on {genome_id}")

            payload = raw.get("model_payload") or {}
            _assert_sequence_free(payload, genome_id)

            receipt = frame_coverage_receipt(raw)
            if receipt is not None:
                frames.append(receipt)
                for seg in receipt.get("segments", []) or []:
                    cid = seg["contig_id"]
                    state = contig_state.setdefault(
                        cid,
                        {
                            "contig_id": cid,
                            "protein_count": 0,
                            "segment_count": 0,
                            "complete": False,
                        },
                    )
                    state["protein_count"] = max(
                        state["protein_count"], int(seg.get("protein_offset_end", 0))
                    )
                    state["segment_count"] += 1
                    state["complete"] = bool(seg.get("complete")) or state["complete"]

                found = self._scan_frame(raw, payload, params, task_id)
                candidates.extend(found)
                LOGGER.info(
                    "  frame %s: %s contigs, %s bytes -> %s candidates",
                    receipt["frame_index"],
                    len(receipt.get("contig_ids", [])),
                    receipt.get("model_payload_bytes"),
                    len(found),
                )

            frame_index += 1
            cursor = raw.get("next_cursor")
            complete = bool(raw.get("complete"))

            if frame_index % checkpoint_interval == 0 or complete:
                self._heartbeat_and_checkpoint(
                    task_id,
                    checkpoint_key,
                    cursor,
                    frames,
                    candidates,
                    contig_state,
                )
            if complete or cursor is None:
                break

        self._finalize(task, params, frames, candidates, contig_state, campaign_id)

    # ------------------------------------------------------------- frame model

    def _scan_frame(
        self,
        raw: dict[str, Any],
        payload: dict[str, Any],
        params: dict[str, Any],
        task_id: str,
    ) -> list[dict[str, Any]]:
        """Step 3 — one bounded model turn per frame."""
        if self.dry_run:
            return []

        run: ModelRun = run_profile(
            provider=self.profile.provider,
            model=self.profile.model,
            reasoning_effort=self.profile.reasoning_effort,
            system_prompt=SCAN_SYSTEM_PROMPT,
            payload_text=_canonical(payload),
            output_schema=SCAN_OUTPUT_SCHEMA,
            timeout=self.model_timeout,
        )

        out: list[dict[str, Any]] = []
        for item in run.record.get("candidates") or []:
            if not isinstance(item, dict):
                continue
            signature = item.get("signature")
            subject_refs = item.get("subject_refs")
            if not isinstance(signature, dict) or not isinstance(subject_refs, dict):
                LOGGER.warning("dropping malformed candidate on frame %s", raw.get("frame_index"))
                continue
            ctype = str(item.get("candidate_type") or "unclassified")
            out.append(
                {
                    "candidate_type": f"{CANDIDATE_TYPE_PREFIX}:{ctype}"[:256],
                    "signature_schema": SIGNATURE_SCHEMA,
                    "signature": signature,
                    "evidence": item.get("evidence") or {},
                    "subject_refs": subject_refs,
                    "reason_codes": list(item.get("reason_codes") or []),
                    "uncertainty": item.get("uncertainty") or {},
                    # Step 8 — provenance travels with the record.
                    "provenance": {
                        "frame_id": raw.get("frame_id"),
                        "frame_index": raw.get("frame_index"),
                        "model_payload_sha256": raw.get("model_payload_sha256"),
                        "packing_contract_hash": raw.get("packing_contract_hash"),
                        "policy_hash": self.policy.policy_hash,
                        "controller_id": self.policy.controller_id,
                        "execution_profile": self.profile_name,
                        "task_id": task_id,
                        **run.provenance(),
                    },
                }
            )
        return out

    # -------------------------------------------------------- fence operations

    def _resume(
        self, task_id: str, checkpoint_key: str
    ) -> tuple[str | None, list[dict], list[dict], dict[str, dict]]:
        """Step 5 (read side) — resume mid-genome across attempts."""
        try:
            ck = self.ops.get_task_checkpoint(task_id, checkpoint_key)
        except requests.HTTPError as exc:
            resp = getattr(exc, "response", None)
            if resp is not None and resp.status_code == 404:
                return None, [], [], {}
            raise
        payload = ck.get("payload") or {}
        frames = payload.get("frames") or []
        LOGGER.info("resuming task %s from frame %s", task_id, len(frames))
        return (
            ck.get("cursor"),
            frames,
            payload.get("candidates") or [],
            payload.get("contig_state") or {},
        )

    def _heartbeat_and_checkpoint(
        self,
        task_id: str,
        checkpoint_key: str,
        cursor: str | None,
        frames: list[dict],
        candidates: list[dict],
        contig_state: dict[str, dict],
    ) -> None:
        try:
            self.ops.heartbeat_task(task_id, lease_seconds=self.lease_seconds)
            self.ops.put_task_checkpoint(
                task_id,
                checkpoint_key,
                cursor=cursor,
                payload={
                    "frames": frames,
                    "candidates": candidates,
                    "contig_state": contig_state,
                },
            )
        except requests.HTTPError as exc:
            if _is_conflict(exc):
                raise LeaseLost(str(exc)) from exc
            raise

    # ---------------------------------------------------------------- finalize

    def _finalize(
        self,
        task: dict[str, Any],
        params: dict[str, Any],
        frames: list[dict],
        candidates: list[dict],
        contig_state: dict[str, dict],
        campaign_id: str | None,
    ) -> None:
        """Steps 6 and 7 — exact coverage, candidates, one disposition, complete."""
        task_id = task["id"]
        unit_id = params["unit_id"]
        genome_id = params["genome_id"]
        dataset_id = params["dataset_id"]

        manifest_path = Path(params["coverage_manifest"])
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest = write_genome_coverage_manifest(
            params, list(contig_state.values()), frames, manifest_path
        )
        coverage_hash = manifest["coverage_sha256"]
        status = manifest.get("coverage_status")
        LOGGER.info(
            "coverage %s: %s/%s contigs, %s/%s proteins, %s frames",
            status,
            manifest.get("observed_n_contigs"),
            manifest.get("expected_n_contigs"),
            manifest.get("observed_n_proteins"),
            manifest.get("expected_n_proteins"),
            manifest.get("model_frame_count"),
        )

        try:
            # Ordering is mandated: every candidate, then exactly one disposition.
            for idx, cand in enumerate(candidates):
                self.ops.create_candidate_occurrence(
                    campaign_id=campaign_id,
                    dataset_id=dataset_id,
                    unit_id=unit_id,
                    genome_id=genome_id,
                    candidate_type=cand["candidate_type"],
                    signature_schema=cand["signature_schema"],
                    signature=cand["signature"],
                    evidence=cand["evidence"],
                    verification=[],
                    subject_refs=cand["subject_refs"],
                    task_id=task_id,
                    reason_codes=cand["reason_codes"],
                    uncertainty=cand["uncertainty"],
                    provenance=cand["provenance"],
                    idempotency_key=_stable_key(
                        task_id, str(idx), _canonical(cand["signature"])
                    ),
                )

            disposition = "candidate" if candidates else "clear"
            if status != "complete":
                disposition = "incomplete"

            self.ops.record_unit_disposition(
                campaign_id=campaign_id,
                unit_id=unit_id,
                dataset_id=dataset_id,
                genome_id=genome_id,
                coverage_hash=coverage_hash,
                candidate_count=len(candidates),
                disposition=disposition,
                evidence_bundle_hash=hashlib.sha256(
                    _canonical([c["signature"] for c in candidates]).encode("utf-8")
                ).hexdigest(),
                task_id=task_id,
                provenance={
                    "policy_hash": self.policy.policy_hash,
                    "controller_id": self.policy.controller_id,
                    "execution_profile": self.profile_name,
                    "agent_id": self.agent_id,
                    "coverage_manifest": str(manifest_path),
                    "frames": len(frames),
                },
                idempotency_key=_stable_key(task_id, "disposition", coverage_hash),
            )

            if disposition == "incomplete":
                # Completion is gated on an active clear/candidate disposition,
                # so an incomplete genome must go back for another attempt
                # rather than being marked done.
                raise RuntimeError(
                    f"coverage incomplete for {genome_id}: {manifest.get('errors')}"
                )

            self.ops.complete_task(task_id)
            LOGGER.info(
                "task %s complete: %s candidates, disposition=%s", task_id, len(candidates), disposition
            )
        except requests.HTTPError as exc:
            if _is_conflict(exc):
                raise LeaseLost(str(exc)) from exc
            raise
