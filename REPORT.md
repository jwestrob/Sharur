## Test Report (current iteration)

Commands run (with `--override-ini addopts=""` to bypass missing pytest-cov plugin):

1. `pytest -q tests/test_tools/test_smoke.py --override-ini addopts=""`  
   - Result: **pass** (2 tests, `datetime.utcnow` deprecation warnings).

2. `pytest -q tests/test_session.py --override-ini addopts=""`  
   - Result: **pass** (3 tests after provenance logging addition, same warnings).

3. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py --override-ini addopts=""`  
   - Result: **pass** (8 tests).  

4. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py tests/test_export.py --override-ini addopts=""`  
   - Result: **pass** (10 tests).  

5. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py tests/test_export.py tests/test_agent_params.py --override-ini addopts=""`  
   - Result: **pass** (11 tests).  

6. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py tests/test_export.py tests/test_agent_params.py --override-ini addopts=""`  
   - Result: **pass** (13 tests after additional export/provenance coverage).  

7. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py tests/test_export.py tests/test_agent_params.py --override-ini addopts=""`  
   - Result: **pass** (14 tests after focus-fallback param extraction).  

8. `pytest -q tests/test_agent_params.py --override-ini addopts=""`  
   - Result: **pass** (focus_hint wiring validation).

9. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py tests/test_export.py tests/test_agent_params.py --override-ini addopts=""`  
   - Result: **pass** (15 tests after working-set param fill).  

10. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py tests/test_export.py tests/test_agent_params.py --override-ini addopts=""`  
    - Result: **pass** (19 tests after scope-fill, export FASTA, multi-tool provenance, exclude_set).  

11. `pytest -q tests/test_session.py tests/test_tools/test_smoke.py tests/test_export.py tests/test_agent_params.py tests/test_ingest_smoke.py --override-ini addopts=""`  
    - Result: **pass** (20 tests including ingest smoke).  

12. `pytest -q tests/test_agent_params.py --override-ini addopts=""`  
    - Result: **pass** (6 tests; expanded DSPy param coverage for anomalies, loci, manage_sets, export).  

13. `pytest -q tests/test_ingest_integration_dataset.py --override-ini addopts=""`  
    - Result: **pass** (integration build against `dummy_dataset` genomes; validates stage07 loading bins/contigs/proteins, annotations, loci, feature_store counts).

14. `pytest -q tests/test_toolchain_on_ingest_dataset.py --override-ini addopts=""`  
    - Result: **pass** (runs core tools against DuckDB built from `dummy_dataset`: find_proteins with PFAM hit, get_context from focus, manage_sets create, export TSV).

Notes:
- Installed `duckdb` locally to satisfy storage backend import during tests.
- Coverage options from pyproject are disabled in these runs; full suite pending dependency setup.
