"""Sharur's shipped KEGG rules and the local build, without KEGG data.

The synthetic build uses invented KOs (K9xxxx) in KEGG's file formats, so it
exercises ``sharur setup-kegg`` end to end with no KEGG content.
"""

import json
import re
from pathlib import Path

import pytest

from sharur.predicates.mappings import kegg_map
from sharur.predicates.mappings.kegg_build import build, load_swissprot_consensus
from sharur.predicates.mappings.kegg_evidence import (
    load_brite_rules,
    load_hyddb_snapshot,
    load_module_rules,
    parse_module_definition,
)
from sharur.predicates.mappings.kegg_map import (
    EC_TO_PREDICATES,
    get_predicates_for_ec,
    get_predicates_for_kegg,
)
from sharur.predicates.mappings.pfam_evidence import (
    SWISSPROT_COVERAGE,
    SWISSPROT_EXCLUDED,
    SWISSPROT_MIN_LOWER_BOUND,
    wilson_lower_bound,
)
from sharur.predicates.vocabulary import ALL_PREDICATES, PREDICATE_BY_ID


DATA = Path(kegg_map.__file__).with_name("data")


def _rows(name):
    return [line.split("\t") for line in (DATA / name).read_text().splitlines()
            if line and not line.startswith("#")]


# --------------------------------------------------------------------------- #
# Shipped rules
# --------------------------------------------------------------------------- #


def test_vocabulary_ids_are_unique():
    ids = [p.predicate_id for p in ALL_PREDICATES]
    assert len(ids) == len(set(ids))


def test_ec_map_names_vocabulary_predicates():
    assert not {p for preds in EC_TO_PREDICATES.values() for p in preds} - set(PREDICATE_BY_ID)


def test_rules_name_vocabulary_predicates():
    named = {p for r in load_brite_rules() for p in r.predicates}
    named |= {p for preds in load_module_rules().values() for p in preds}
    named |= {p for _, preds, *_ in _rows("kegg_predicate_proposals.tsv") for p in preds.split(",")}
    named |= {pred for _, pred, *_ in _rows("kegg_swissprot_consensus.tsv")}
    assert not named - set(PREDICATE_BY_ID)


def test_repository_ships_no_kegg_text():
    # KEGG names, BRITE placements and module definitions are built locally by `sharur setup-kegg`.
    assert not (DATA / "kegg_predicates.tsv").exists()
    assert not (DATA / "kegg_evidence_snapshot.tsv").exists()
    assert all(len(row) == 2 and re.fullmatch(r"M\d{5}", row[0]) for row in _rows("kegg_module_predicates.tsv"))
    symbol = re.compile(r"[\w.\-/+()']*")
    assert all(len(row) == 3 and symbol.fullmatch(row[2]) for row in _rows("kegg_predicate_proposals.tsv"))


def test_shipped_consensus_meets_the_thresholds():
    failures = []
    for ko, pred, k, n in _rows("kegg_swissprot_consensus.tsv"):
        k, n = int(k), int(n)
        ok = pred not in SWISSPROT_EXCLUDED and n >= SWISSPROT_COVERAGE[0] and k >= SWISSPROT_COVERAGE[1] * n \
            and wilson_lower_bound(k, n) >= SWISSPROT_MIN_LOWER_BOUND
        if not ok:
            failures.append((ko, pred, k, n))
    assert not failures, failures[:10]


def test_hyddb_snapshot_counts_are_consistent():
    for ko, labels in load_hyddb_snapshot().items():
        assert all(n > 0 for n in labels.values()), ko


def test_module_definition_parser():
    parsed = parse_module_definition("K00330+(K00331+K00332,K13380)+K00334 (K00844,K12407) K01803-(K1,K00002)")
    assert all(parsed[k] for k in ("K13380", "K00334", "K00002"))
    assert not any(parsed[k] for k in ("K00844", "K12407"))


def test_runtime_without_a_map_uses_kegg_stated_ec_only(monkeypatch):
    monkeypatch.setattr(kegg_map, "KEGG_TO_PREDICATES", {})
    assert get_predicates_for_kegg("K99999", "hydrogenase membrane protein") == []
    assert set(get_predicates_for_kegg("K99999", "x [EC:2.7.1.1]")) == set(get_predicates_for_ec("2.7.1.1"))


@pytest.mark.parametrize(("ec", "present", "absent"), [
    ("2.1.2.1", {"transferase"}, {"methyltransferase"}),
    ("2.1.1.37", {"methyltransferase", "dna_methylase"}, set()),
    ("2.1.1.45", {"methyltransferase"}, {"dna_methylase"}),
    ("3.6.5.2", {"gtp_binding"}, {"atpase"}),
    ("3.6.1.1", {"hydrolase"}, {"atpase", "atp_binding"}),
    ("3.1.21.4", {"nuclease", "endonuclease"}, {"esterase"}),
    ("3.1.-.-", {"hydrolase"}, {"esterase"}),
    ("4.1.3.1", {"lyase"}, {"decarboxylase"}),
    ("7.6.2.1", {"atp_binding"}, {"gtp_binding"}),
    ("2.7.8.5", {"transferase"}, {"kinase"}),
    ("1.4.1.3", {"oxidoreductase"}, {"aminotransferase"}),
    ("2.6.1.1", {"aminotransferase", "plp_binding"}, set()),
    ("1.2.1.12", {"nad_binding"}, set()),                          # ENZYME 1.2.1: NAD(P)+ acceptor
    ("1.3.8.7", {"flavin_binding"}, {"nad_binding"}),              # ENZYME 1.3.8: flavin acceptor
    ("1.14.12.10", {"dioxygenase", "nad_binding"}, {"monooxygenase"}),
    ("1.14.13.1", {"monooxygenase", "nad_binding"}, {"dioxygenase"}),
    ("1.14.11.2", {"dioxygenase"}, {"monooxygenase"}),
    ("1.14.19.1", {"oxygenase"}, {"monooxygenase", "dioxygenase"}),
])
def test_ec_classes(ec, present, absent):
    preds = set(get_predicates_for_ec(ec))
    assert present <= preds
    assert not absent & preds


# --------------------------------------------------------------------------- #
# Synthetic end-to-end build
# --------------------------------------------------------------------------- #


def _synthetic_inputs(root: Path) -> tuple[Path, Path]:
    inputs = root / "inputs"
    (inputs / "brite").mkdir(parents=True)
    (inputs / "ko_list.tsv").write_text(
        "ko:K90001\tfooA; foo kinase [EC:2.7.1.1]\n"
        "ko:K90002\tbarR; two-component system, OmpR family, bar response regulator BarR\n"
        "ko:K90003\tbazA; baz complex subunit A\n"
        "ko:K90004\tbazB; baz complex subunit B\n"
        "ko:K90006\tquxA; unplaced protein\n"
    )
    (inputs / "kegg_info.txt").write_text("kegg\tKEGG\n\tbrite 1 2026/01/02\n\tmodule 1 2026/01/03\n\tko 5 2026/01/01\n")
    for hierarchy in {r.hierarchy for r in load_brite_rules()}:
        tree = {"name": hierarchy, "children": []}
        if hierarchy == "ko02022":
            tree["children"] = [{"name": "OmpR family", "children": [
                {"name": "BarS-BarR (bar sensing)", "children": [{"name": "K90002  barR; bar response regulator"}]}]}]
        (inputs / "brite" / f"{hierarchy}.json").write_text(json.dumps(tree))
    entries = [f"ENTRY       {m}            Pathway   Module\nDEFINITION  K99998\n///\n"
               for m in load_module_rules() if m != "M00144"]
    entries.append("ENTRY       M00144            Pathway   Module\nDEFINITION  K90003+K90004\n///\n")
    (inputs / "modules.txt").write_text("".join(entries))
    kofam = root / "kofam_ko_list"
    kofam.write_text("knum\tthreshold\tscore_type\tdefinition\nK90005\t100\tfull\tretired thing [EC:1.1.1.1]\n")
    return inputs, kofam


def test_synthetic_build_end_to_end(tmp_path, monkeypatch):
    inputs, kofam = _synthetic_inputs(tmp_path)
    out = tmp_path / "kegg"
    provenance = build(inputs, out, kofam_ko_list=kofam, log=lambda *_: None)

    predicates, evidence = kegg_map._load_kegg_mapping_file(out / "kegg_predicates.tsv")
    assert {"kinase", "transferase"} <= set(predicates["K90001"])
    assert evidence["K90001"]["kinase"] == "ec:2.7.1.1"
    assert {"two_component", "response_regulator", "ompr_family"} <= set(predicates["K90002"])
    assert evidence["K90002"]["ompr_family"].startswith("brite:ko02022 OmpR family")
    assert evidence["K90003"]["complex_subunit"] == "module:M00144"
    assert evidence["K90004"]["electron_transport"] == "module:M00144"
    assert {"dehydrogenase", "nad_binding"} <= set(predicates["K90005"])  # KOfam-only KO
    assert "K90006" not in predicates

    assert (out / "kegg_predicates.tsv").read_text().startswith("# Built locally by sharur setup-kegg")
    assert provenance["kegg_release"] == "ko 2026/01/01; brite 2026/01/02; module 2026/01/03"
    assert provenance["map_sha256"] == json.loads((out / "provenance.json").read_text())["map_sha256"]
    assert set(provenance["rules"]) >= {"kegg_brite_predicates.tsv", "kegg_swissprot_consensus.tsv"}

    monkeypatch.setenv(kegg_map.KEGG_DIR_ENV, str(out))
    assert kegg_map.kegg_dir() == out


def test_build_refuses_rules_for_missing_modules(tmp_path):
    inputs, kofam = _synthetic_inputs(tmp_path)
    (inputs / "modules.txt").write_text("ENTRY       M00144            Pathway   Module\nDEFINITION  K90003\n///\n")
    with pytest.raises(ValueError, match="absent from this KEGG release"):
        build(inputs, tmp_path / "kegg", kofam_ko_list=kofam, log=lambda *_: None)


def test_consensus_table_loads_as_evidence():
    table = load_swissprot_consensus()
    ko, pred, k, n = _rows("kegg_swissprot_consensus.tsv")[0]
    assert table[ko][pred] == f"swissprot:{k}/{n}"
