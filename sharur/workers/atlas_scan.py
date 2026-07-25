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
from sharur.workers.model_cli import (
    ModelError,
    ModelRateLimited,
    ModelRun,
    ModelTransient,
    run_profile,
)

LOGGER = logging.getLogger("sharur.workers.atlas_scan")

TASK_TYPE = "atlas_genome_read"
CANDIDATE_TYPE_PREFIX = "atlas-scan"
SIGNATURE_SCHEMA = "atlas-scan-signature/1.0"

# What the model is allowed to emit.
#
# OpenAI strict structured outputs requires every object to set
# `additionalProperties: false` AND to list every property in `required` —
# free-form objects are rejected outright with `invalid_json_schema`. The
# signature/evidence/subject_refs fields are inherently open-ended (the model
# chooses the typed fields that describe a finding), so they travel as JSON
# **strings** and are parsed worker-side. That keeps a hard guarantee on the
# outer shape while leaving the inner content free.
SCAN_OUTPUT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["candidates", "notes"],
    "properties": {
        "candidates": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": [
                    "candidate_type",
                    "signature",
                    "evidence_json",
                    "subject_refs_json",
                    "reason_codes",
                    "confidence",
                ],
                "properties": {
                    "candidate_type": {"type": "string"},
                    # The reduction key is the one field that must be
                    # canonical, so it is a fixed structure rather than
                    # free-form JSON. Left open, the model invents a different
                    # schema per finding (`accessions` vs `caller_accessions`
                    # vs `defense_calls` vs `components`), and the same system
                    # in two genomes hashes differently — observed at 6
                    # competing schemas across 37 signatures in the first run.
                    "signature": {
                        "type": "object",
                        "additionalProperties": False,
                        "required": ["accessions", "n_genes", "system", "subtype"],
                        "properties": {
                            "accessions": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "Every PFAM/KEGG/caller accession in the locus. Order is normalised worker-side.",
                            },
                            "n_genes": {
                                "type": "integer",
                                "description": "Number of genes in the locus",
                            },
                            "system": {
                                "type": "string",
                                "description": (
                                    "COARSE module type — the grouping key. Use the broad "
                                    "category, e.g. glycan-locus, defense-island, "
                                    "secretion-system, bgc. Never compound it with detail."
                                ),
                            },
                            "subtype": {
                                "type": "string",
                                "description": (
                                    "OPTIONAL refinement of `system`, e.g. pul, capsule-eps, "
                                    "t2ss, nrps, cbass-cargo. Empty string when the evidence "
                                    "does not resolve it — that is a normal, useful answer."
                                ),
                            },
                        },
                    },
                    "evidence_json": {
                        "type": "string",
                        "description": "JSON object holding reasoning and supporting observations",
                    },
                    "subject_refs_json": {
                        "type": "string",
                        "description": 'JSON object, e.g. {"protein_ids": [...], "contig_ids": [...]}',
                    },
                    "reason_codes": {"type": "array", "items": {"type": "string"}},
                    "confidence": {
                        "type": "string",
                        "enum": ["low", "medium", "high"],
                    },
                },
            },
        },
        "notes": {"type": "string"},
    },
}


def _parse_json_field(raw: Any, field: str) -> dict[str, Any] | None:
    """Parse one of the JSON-string carrier fields into a dict."""
    if isinstance(raw, dict):
        return raw  # tolerated: the Anthropic path is not schema-constrained
    if not isinstance(raw, str) or not raw.strip():
        return None
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        LOGGER.warning("candidate %s was not valid JSON; dropping candidate", field)
        return None
    return parsed if isinstance(parsed, dict) else None

# Fields that locate a finding within one genome. They must never enter the
# signature: the signature is the exact reduction key, so any genome-specific
# value makes identical biology in two genomes hash differently and defeats
# cross-genome reduction entirely.
#
# Observed in the first pilot run: `contig_id` appeared in 363 of 400
# signatures and the collapse ratio was 1.00x — 398 of 399 clusters were
# singletons. Prompting alone is not enough (models invent key names like
# `gene_indexes`, `gene_index_start`, `locus_end`), so this is enforced by
# pattern rather than by an exact blacklist, and locators are moved into
# subject_refs where they belong.
_LOCATOR_PATTERNS = (
    "contig",
    "protein_id",
    "coordinate",
    "position",
    "offset",
    "locus_start",
    "locus_end",
    "gene_index",
    "gene_indices",
    "gene_indexes",
    "start",
    "end",
    "span",
    "range",
)


def _is_locator(key: str) -> bool:
    k = key.lower()
    return any(pat in k for pat in _LOCATOR_PATTERNS)


# System labels drift ("surface-polysaccharide-biosynthesis" /
# "surface-polysaccharide" / "surface-polysaccharide-locus" are one thing), so
# the name is normalised before it becomes part of the key.
# Coarse module types. These are the grouping keys a human triages — "there are
# ten glycan loci" — so the list is deliberately short. Refinement goes in
# `subtype`, which never enters the reduction key.
_BASE_MODULE_TYPES = frozenset(
    {
        "glycan-locus",
        "defense-island",
        "secretion-system",
        "bgc",
        "mobile-element",
        "prophage",
        "crispr-cas",
        "hydrogenase-system",
        "co-dehydrogenase",
        "formate-dehydrogenase",
        "carbon-fixation",
        "cofactor-biosynthesis",
        "storage-granule",
        "microcompartment",
        "cell-wall-modification",
        "surface-adhesion",
        "motility",
        "natural-competence",
        "stress-response",
    }
)

_SYSTEM_NOISE = (
    "-locus", "_locus", " locus",
    "-biosynthesis", "_biosynthesis", " biosynthesis",
    "-cluster", "_cluster", " cluster",
    "-system", "_system", " system",
    "-operon", "_operon", " operon",
)


def _normalise_system(name: str) -> str:
    """Fold label drift, without damaging a canonical base type.

    Suffix stripping exists to collapse drift like
    "surface-polysaccharide-biosynthesis" onto "surface-polysaccharide". But
    several base types legitimately end in one of those tokens (glycan-locus,
    defense-island, secretion-system), so an exact base-type match short
    circuits the stripping — otherwise `glycan-locus` degrades to `glycan` and
    stops matching the vocabulary the demotion logic checks against.
    """
    out = (name or "").strip().lower().replace("_", "-").replace(" ", "-")
    if out in _BASE_MODULE_TYPES:
        return out
    changed = True
    while changed:
        changed = False
        for suffix in _SYSTEM_NOISE:
            token = suffix.strip().replace("_", "-").replace(" ", "-")
            if token and out.endswith(token) and len(out) > len(token):
                candidate = out[: -len(token)].rstrip("-")
                if candidate in _BASE_MODULE_TYPES:
                    return candidate
                out = candidate
                changed = True
    return out


def _canonical_signature(
    sig: dict[str, Any], candidate_type: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Split a model signature into (reduction key, reduction features).

    The key must be the *coarse identity* of a finding, not its exact gene
    complement. Two genuine Gabija systems differ by a domain or two, so an
    exact match on the full accession set never collides: measured on live
    output, sharing fell from 16% of clusters at five accessions to 2% at ten,
    giving a 1.12x collapse ratio with 94% singletons.

    The named system, by contrast, repeats strongly (Gabija x23, Mokosh type I
    x20, Hachiman x12 across nineteen genomes), so it is the right granularity.
    The accession set is preserved in `reduction_features`, where the reducer
    and reviewers can see the per-occurrence variation without it fragmenting
    the key.

    Occurrences with no named system fall back to the accession set, since a
    coarse key of "unnamed" would wrongly merge unrelated findings.
    """
    accessions = sig.get("accessions") or []
    if not isinstance(accessions, list):
        accessions = [accessions]
    accessions = sorted({str(a) for a in accessions})
    n_genes = int(sig.get("n_genes") or 0)
    system = _normalise_system(str(sig.get("system") or ""))
    subtype = _normalise_system(str(sig.get("subtype") or ""))
    cls = candidate_type.replace(CANDIDATE_TYPE_PREFIX + ":", "")

    # A compound label fragments the very groups it should form: emitted
    # free-form, five defence findings arrived as five singleton labels
    # (defense-island-cbass-pmt-cargo, -mokosh-signaling-array, ...) instead of
    # one group of eight. Detail after the base type is demoted to `subtype`.
    if system.count("-") >= 2 and not system.startswith("novel"):
        head, _, tail = system.partition("-")
        second, _, rest = tail.partition("-")
        base = f"{head}-{second}"
        if base in _BASE_MODULE_TYPES:
            system, subtype = base, subtype or f"{rest}"
        elif head in _BASE_MODULE_TYPES:
            system, subtype = head, subtype or tail

    if system:
        key = {"system": system, "class": cls}
    else:
        key = {"system": "", "class": cls, "accessions": accessions}
    features = {"accessions": accessions, "n_genes": n_genes, "subtype": subtype}
    return key, features


def _split_signature(
    signature: dict[str, Any], genome_id: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Partition a model signature into (reduction key, locators).

    Anything naming a position, or whose value embeds the genome id, is a
    locator. What remains describes *what the finding is* — accession sets,
    counts, typed classes — and reduces across genomes.
    """
    invariant: dict[str, Any] = {}
    locators: dict[str, Any] = {}
    for key, value in signature.items():
        if _is_locator(key) or (genome_id and genome_id in _canonical(value)):
            locators[key] = value
        else:
            invariant[key] = value
    return invariant, locators


SCAN_SYSTEM_PROMPT = """\
You are an Atlas scanner reading one bounded frame of one microbial genome.

You will receive a JSON packet containing contigs from exactly one genome, each
with its proteins. Every protein carries: protein_id, gene_index, start, end,
strand, length_aa, gc_content, observed_annotations, predicates, named_calls,
and loci. There are no biological sequences and you must never ask for them.

Your job is to report NOTABLE architecture as typed candidate occurrences. Do
not summarize the genome.

YOUR OUTPUT IS AN INDEX, NOT AN ANALYSIS. Findings are grouped downstream by
the `system` label you assign, and a human triages the groups — "there are ten
of this BGC class" — then decides which deserve a closer look. So:

MARK LIBERALLY. A marked-but-ordinary locus costs one glance. A missed one is
lost for the whole campaign. Do not suppress a borderline call.

LABEL PRECISELY. The `system` label is the single most important thing you
emit: it is both the grouping key and the unit of triage. A vague label buries
distinct biology in one undifferentiated pile.

Assign a COARSE `system` and, when the evidence supports it, a `subtype`.

`system` is the grouping key. Keep it broad and never compound it with detail:

  glycan-locus, defense-island, secretion-system, bgc, mobile-element,
  prophage, crispr-cas, hydrogenase-system, co-dehydrogenase,
  formate-dehydrogenase, carbon-fixation, cofactor-biosynthesis,
  storage-granule, microcompartment, cell-wall-modification,
  surface-adhesion, motility, natural-competence, stress-response

`subtype` refines it, and is OPTIONAL:

  glycan-locus     -> pul | capsule-eps | o-antigen | protein-glycosylation
  secretion-system -> t2ss | t3ss | t4ss | t6ss | t7ss
  bgc              -> nrps | pks | ripp | terpene | hybrid
  defense-island   -> whatever the cargo or architecture is, e.g. cbass-cargo
  carbon-fixation  -> cbb | wood-ljungdahl | 3hp

**Leave `subtype` empty when the evidence does not resolve it.** A glycan locus
you cannot assign to PUL versus capsule is still a real, reportable glycan
locus — record it as `system: glycan-locus`, `subtype: ""`, and say in
`evidence` which discriminating feature was missing (a SusC/SusD-like pair, a
Wzy/Wzx export system, and so on). An empty subtype is a normal and useful
answer; it groups cleanly and can be revisited. Do NOT guess a subtype to avoid
leaving it blank, and do NOT invent a compound `system` to smuggle detail into
the key.

Use `novel-<something>` as the `system` only when the locus fits no category
above at all.

Naming a module is a CLAIM, and a claim can be wrong. That is intended. A
finding that cannot be wrong ("these genes are glycosyltransferases, therefore
a glycosyltransferase locus") restates the input and says nothing. State what
the module IS and what it DOES.

DO NOT RE-INVENTORY WHAT A PURPOSE-BUILT CALLER ALREADY PRODUCED. Defence and
antiviral systems are already called by DefenseFinder, with gene membership,
and they reach you in `named_calls`. Never emit a candidate whose content is
"this genome has a Gabija/Mokosh/Hachiman/RM/CBASS system" — that is known, and
better established than anything inferable from adjacent PFAM domains. Report a
defence-associated locus ONLY when the CONTEXT is the finding — the system sits
inside a mobile element, carries unusual cargo, or is adjacent to
uncharacterised genes suggesting an extension — and say so explicitly.

DO NOT report a single annotated gene with no informative neighbourhood, or
core machinery whose only claim is that it exists (ribosome, ATP synthase,
Complex I, chaperones, the dcw/cell-division cluster, lone transporters or
two-component systems).

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

Return one JSON object with a "candidates" array and a "notes" string. Three
fields carry a JSON **object encoded as a string** — write valid JSON inside
the string:

  candidate_type      short kebab-case class, e.g. "co-located-pathway",
                      "novel-gene-cluster", "multidomain-architecture",
                      "mobile-element-cargo", "defense-locus"
  signature_json      stringified JSON object describing WHAT THE FINDING IS,
                      never where it is. Accession sets, counts, typed classes
                      only. This is the exact reduction key: the SAME finding in
                      a different genome must produce a BYTE-IDENTICAL
                      signature, so it must contain no contig id, no
                      coordinates, no gene indices, no protein ids and no
                      genome name. Put those in subject_refs_json instead.
                      e.g. "{\\"accessions\\":[\\"PRK\\",\\"RuBisCO_large\\",\\"RuBisCO_small\\"],\\"n_genes\\":3}"
  evidence_json       stringified JSON object holding your reasoning, the
                      supporting observations, and any UNVERIFIED hypotheses
  subject_refs_json   stringified JSON object, e.g.
                      "{\\"protein_ids\\":[\\"p1\\"],\\"contig_ids\\":[\\"c1\\"]}"
  reason_codes        array of short tags (may be empty)
  confidence          one of "low", "medium", "high"

All six fields are required on every candidate. "notes" may be an empty string.

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
        transient_retries: int = 4,
        sweep_failed: bool = True,
        max_sweeps: int = 20,
        dry_run: bool = False,
    ) -> None:
        self.agent_id = agent_id
        self.profile_name = profile
        self.campaign_id = campaign_id
        self.lease_seconds = lease_seconds
        self.model_timeout = model_timeout
        self.max_frames = max_frames
        self.transient_retries = transient_retries
        self.sweep_failed = sweep_failed
        self.max_sweeps = max_sweeps
        self._sweeps = 0
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
                # Queue looks drained. Before idling, go back for tasks that
                # exhausted their attempts — on this cluster that is usually
                # transport (intermittent DNS), not a defect in the genome.
                if self._sweep_failed_tasks():
                    continue
                time.sleep(idle_sleep)
                continue

            try:
                self.run_task(task)
                done += 1
                rate_limit_backoff = 0.0
            except LeaseLost:
                LOGGER.warning("lease lost on task %s; abandoning", task.get("id"))
            except ModelTransient as exc:
                # Exhausted in-process retries: the network is genuinely down.
                # Release for retry and wait rather than spinning on claims.
                LOGGER.warning(
                    "transport still failing after %s retries; releasing %s and sleeping 120s: %s",
                    self.transient_retries,
                    task.get("id"),
                    str(exc)[:200],
                )
                self._release(task, f"transient transport: {exc}", retry_delay=120)
                time.sleep(120)
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

    def _sweep_failed_tasks(self) -> bool:
        """Requeue attempt-exhausted tasks. Returns True if any were reset.

        Only failures that read as transport are reset, so a genuinely broken
        genome is not retried forever. Bounded by `max_sweeps` as a backstop.

        Several workers may sweep at once; the underlying UPDATE is atomic and
        filtered on `status = 'failed'`, so the loser of a race simply resets
        nothing. A short jittered pause keeps them from thundering.
        """
        if not self.sweep_failed or self._sweeps >= self.max_sweeps:
            return False
        time.sleep(random.uniform(0, 3.0))
        try:
            result = self.ops.reset_failed_tasks(
                campaign_id=self.campaign_id,
                only_transient=True,
            )
        except requests.HTTPError as exc:
            LOGGER.warning("failed-task sweep rejected: %s", exc)
            return False
        reset = result.get("reset") or []
        skipped = result.get("skipped") or []
        if not reset:
            if skipped:
                LOGGER.info(
                    "sweep: %s failed task(s) left alone (non-transient errors)",
                    len(skipped),
                )
            return False
        self._sweeps += 1
        LOGGER.info(
            "sweep %s/%s: requeued %s attempt-exhausted task(s); %s left alone",
            self._sweeps,
            self.max_sweeps,
            len(reset),
            len(skipped),
        )
        return True

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

        # Transient transport failures (DNS, socket resets) are retried in
        # place. Releasing the task instead would consume an attempt, and at
        # the default max_attempts=5 a handful of DNS blips would permanently
        # kill a genome that has nothing wrong with it.
        last: Exception | None = None
        run: ModelRun | None = None
        for attempt in range(self.transient_retries + 1):
            try:
                run = run_profile(
                    provider=self.profile.provider,
                    model=self.profile.model,
                    reasoning_effort=self.profile.reasoning_effort,
                    system_prompt=SCAN_SYSTEM_PROMPT,
                    payload_text=_canonical(payload),
                    output_schema=SCAN_OUTPUT_SCHEMA,
                    timeout=self.model_timeout,
                )
                break
            except ModelTransient as exc:
                last = exc
                if attempt >= self.transient_retries:
                    break
                delay = min(2.0 * (2**attempt), 60.0) + random.uniform(0, 2.0)
                LOGGER.warning(
                    "transient transport failure on frame %s (attempt %s/%s): %s; retrying in %.0fs",
                    raw.get("frame_index"),
                    attempt + 1,
                    self.transient_retries,
                    str(exc)[:200],
                    delay,
                )
                self._keepalive(task_id)
                time.sleep(delay)
        if run is None:
            raise ModelTransient(
                f"transport failed {self.transient_retries + 1}x on frame "
                f"{raw.get('frame_index')}: {last}"
            )

        out: list[dict[str, Any]] = []
        for item in run.record.get("candidates") or []:
            if not isinstance(item, dict):
                continue
            signature = item.get("signature")
            if not isinstance(signature, dict):
                signature = _parse_json_field(item.get("signature_json"), "signature")
            subject_refs = _parse_json_field(
                item.get("subject_refs_json", item.get("subject_refs")), "subject_refs"
            )
            evidence = _parse_json_field(
                item.get("evidence_json", item.get("evidence")), "evidence"
            ) or {}
            if signature is None or subject_refs is None:
                LOGGER.warning("dropping malformed candidate on frame %s", raw.get("frame_index"))
                continue
            ctype = str(item.get("candidate_type") or "unclassified")
            genome_id = raw.get("bin_id") or ""
            signature, locators = _split_signature(signature, genome_id)
            signature, reduction_features = _canonical_signature(signature, ctype)
            if not signature:
                LOGGER.warning(
                    "candidate on frame %s had no genome-invariant signature fields; dropping",
                    raw.get("frame_index"),
                )
                continue
            if locators:
                subject_refs = {**subject_refs, "locators": locators}
            out.append(
                {
                    "candidate_type": f"{CANDIDATE_TYPE_PREFIX}:{ctype}"[:256],
                    "signature_schema": SIGNATURE_SCHEMA,
                    "signature": signature,
                    "evidence": evidence,
                    "subject_refs": subject_refs,
                    "reason_codes": list(item.get("reason_codes") or []),
                    "uncertainty": (
                        {"confidence": item["confidence"]}
                        if item.get("confidence")
                        else {}
                    ),
                    # Step 8 — provenance travels with the record.
                    "reduction_features": reduction_features,
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

    def _keepalive(self, task_id: str) -> None:
        """Refresh the lease during a retry backoff.

        Without this a long transport outage would let the lease expire mid
        genome, the coordinator would reclaim the task, and our eventual write
        would 409 against a dead fence.
        """
        try:
            self.ops.heartbeat_task(task_id, lease_seconds=self.lease_seconds)
        except requests.HTTPError as exc:
            if _is_conflict(exc):
                raise LeaseLost(str(exc)) from exc

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
                    reduction_features=cand.get("reduction_features") or {},
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
