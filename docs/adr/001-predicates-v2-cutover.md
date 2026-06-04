# ADR-001: Predicates V2 — Cutover, Shadow, or Archive?

**Status:** Proposed — awaiting decision
**Date:** 2026-06-03
**Decision needed by:** End of next ingest wave (otherwise the question gets stale)

## Context

The repo carries two parallel predicate systems:

- **V1** (`sharur/predicates/`) — flat boolean predicates. Used by every skill.
- **V2** (`sharur/predicates_v2/`) — typed semantic atoms with facets, relations, evidence strength, composite YAML rules, and an unmapped-accession review queue.

V2 was built to run in shadow mode alongside V1. The infrastructure is complete:
`generate_and_persist_v2`, `shadow_diff`, `build_review_queue`, and the compat
layer `semantic_state_to_predicates` that maps V2 back to V1's flat predicate
shape so downstream callers don't need to change.

What does NOT exist:

- A skill that calls V2 directly.
- A coordinator dispatch that uses V2 atoms instead of V1 predicates.
- Recent shadow-diff runs comparing V1 and V2 outputs on a real dataset.
- A target dataset that has both V1 and V2 tables populated.

V2 has not been touched since the schema was built. It is **complete code, untested in live use, and unreferenced from any agent skill.**

## What V2 would actually fix

The biggest verification gap (per the recent system review) is annotation
laundering: PFAM `SAVED_N` gets reported as "Mokosh activity"; PF13437
`MtrB_PioB` gets reported as "metal-reducing protein"; KOs at the superfamily
boundary get cited as the iconic member function.

V1 cannot catch this — boolean predicates compress all evidence levels into
"protein has predicate X = True." V2's facet model directly addresses it:

- `(activity, hydrogenase, implies)` from a HydDB Mbh hit is strong.
- `(architecture, fe_only_cluster, supports)` from a PFAM domain is weaker.
- `(activity, fefe_hydrogenase, supports)` from PFAM alone — composite rule
  refuses to assert activity without a corroborating atom.

The composite YAML in `config/predicates_v2/composites.yaml` is the leverage
point — declarative rules that say "do not assert X without evidence Y." This
is exactly the layer that closes the iconic-member trap.

## Options

### Option A: Two-week shadow comparison, then decide

Pick two datasets with mature V1 annotations (`coronamine_boiler_100nm` and
`spicy_lams`). Run `generate_and_persist_v2` on both. Run `shadow_diff` and
measure:

- How many V1 predicates does V2 not generate? (V2 too strict?)
- How many V2 atoms have no V1 equivalent? (V2 captures more?)
- For high-novelty findings already published in coronamine — would V2 have
  generated the same predicate set?
- What does `build_review_queue` surface? Are the unmapped accessions the
  ones we know cause confusion (SAVED_N, MnhA paralogs of MbhD, etc.)?

Then decide cutover vs continued shadow with hard evidence.

**Cost:** 2 weeks of patience, modest compute on biotite. No agent rewrites
during the window — skills keep using V1.

**Risk:** If V2 has structural bugs, we find out now; if not, we move forward
with confidence.

### Option B: Commit to a phased cutover

Pick one new ingest dataset (`srvp_bacteria_pb` is the obvious candidate since
it's actively ingesting) as the V2 testbed. Add a `generate_v2: bool = True`
default to `sharur-ingest` so V2 atoms get populated at build time. Update
exactly one skill — `defense.md` is the right candidate (the false-positive
rate problem is documented) — to call V2 atoms when present.

Run V1 and V2 in parallel for 4 weeks. If the V2-using defense skill produces
materially fewer false positives, expand to prophage / hydrogenase.

**Cost:** Schema cost in DuckDB (small), one skill rewrite (~half-day), 4
weeks of evaluation.

**Risk:** Committing to a path before measuring. If V2 misses something V1
catches, the new defense skill could regress.

### Option C: Archive V2

Move `sharur/predicates_v2/` and `config/predicates_v2/` to `_attic/`. Write
an ADR documenting the design and why it was set aside (probably: V1 works,
labor cost of skill rewrites is high, the iconic-member trap is being
addressed prose-side via `_validation_protocols.md`).

**Cost:** Small. Keeps the codebase focused.

**Risk:** The actual problem V2 was designed to fix is unsolved, and we'll
hit it again next analysis.

## Recommendation

**Option A.** Reason: V2 was built specifically for the problem that hurts us
most, but we've never measured whether it actually solves it on real data.
Two weeks of shadow-comparison is the cheapest way to find out.

If the shadow diff shows V2 catching real V1 errors → go to Option B with
evidence. If it shows V2 just reproducing V1's outputs → archive cleanly. If
it shows V2 too strict (refuses to fire on cases V1 caught correctly) →
tune the composites, repeat.

The worst outcome is no decision — V2 stays half-integrated forever, V1 keeps
producing laundering errors, the codebase carries dead weight.

## Concrete next steps if Option A is approved

1. Run `b.generate_v2()` on `data/coronamine_boiler_100nm/sharur.duckdb`.
2. Run `b.generate_v2()` on `data/spicy_lams/sharur.duckdb`.
3. Run `shadow_diff` for both. Persist the diff to
   `data/<dataset>/predicates_v2_diff.json`.
4. Run `build_review_queue` on each. Manually triage the top ~50 unmapped
   accessions per dataset.
5. Take 5 high-novelty findings from each dataset's `exploration/findings.jsonl`.
   For each finding's predicate-derived claim, ask: does V2 generate the same
   atom set? Does the composite rule fire?
6. Write `docs/adr/002-predicates-v2-decision.md` documenting the result.

Time budget: half-day for steps 1-3 (compute), half-day for 4-5 (triage),
half-day for 6 (write-up). Two days of focused work spread over two weeks
keeps it from blocking other work.

## Sign-off

- [ ] Approved — proceed with Option A
- [ ] Approved — proceed with Option B (specify testbed dataset)
- [ ] Approved — proceed with Option C (archive)
- [ ] Deferred — revisit at end of next ingest wave

_Decision: ___________________________________  _Date: __________________
