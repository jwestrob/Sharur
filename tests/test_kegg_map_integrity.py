"""KEGG/EC predicate maps agree with KEGG and with the predicate vocabulary.

`data/kegg_ko_names.tsv` snapshots KEGG's name for every mapped KO (the KOfam
ko_list definition for KOs since retired from KEGG REST). Each map
entry's comment starts with KEGG's symbol for that KO, so a mistyped or shifted
KO number fails here instead of silently labelling a different protein.
"""

import re
from pathlib import Path

import pytest

from sharur.predicates.mappings import kegg_map
from sharur.predicates.mappings.kegg_map import (
    EC_TO_PREDICATES,
    KEGG_TO_PREDICATES,
    get_predicates_for_ec,
)
from sharur.predicates.vocabulary import PREDICATE_BY_ID


SNAPSHOT = Path(kegg_map.__file__).with_name("data") / "kegg_ko_names.tsv"
SOURCE = Path(kegg_map.__file__).read_text()
ENTRY = re.compile(r'^\s*"(K\d{5})":\s*\[[^\]]*\],?\s*#\s*(.*)$', re.M)


def _snapshot() -> dict[str, tuple[set[str], str]]:
    rows = {}
    for line in SNAPSHOT.read_text().splitlines():
        if line.startswith("#"):
            continue
        ko, symbols, definition = line.split("\t")
        rows[ko] = ({s.strip().lower() for s in symbols.split(",")}, definition)
    return rows


def test_every_mapped_ko_is_in_the_kegg_snapshot():
    assert set(KEGG_TO_PREDICATES) <= set(_snapshot())


def test_every_comment_starts_with_the_kegg_symbol():
    snapshot = _snapshot()
    comments = dict(ENTRY.findall(SOURCE))
    assert set(comments) == set(KEGG_TO_PREDICATES)
    mismatched = {
        ko: comment for ko, comment in comments.items()
        if comment.split(";")[0].strip().lower() not in snapshot[ko][0]
    }
    assert not mismatched


@pytest.mark.parametrize("table", [KEGG_TO_PREDICATES, EC_TO_PREDICATES])
def test_every_mapped_predicate_is_in_the_vocabulary(table):
    unknown = {p for preds in table.values() for p in preds if p not in PREDICATE_BY_ID}
    assert not unknown


@pytest.mark.parametrize(("ko", "present", "absent"), [
    ("K03529", {"chromosome_partitioning"}, {"crispr_associated", "defense_system"}),  # smc
    ("K03595", {"gtp_binding"}, {"chromosome_partitioning"}),                           # era
    ("K01356", {"repressor"}, {"integrase", "mobile_element"}),                         # lexA
    ("K15830", {"formate_coupled"}, {"ech_hydrogenase"}),                               # hycE
    ("K14087", {"ech_hydrogenase"}, set()),                                             # echB
    ("K14126", {"nife_group3"}, {"mbh_hydrogenase", "nife_group4"}),                    # mvhA
    ("K18008", {"nife_hydrogenase"}, {"fefe_hydrogenase"}),                             # [NiFe] hydA
    ("K00436", {"nad_coupled"}, {"uptake_hydrogenase"}),                                # hoxH
    ("K02588", {"nitrogenase"}, set()),                                                 # nifH
    ("K02040", {"phosphate_transporter"}, {"sulfate_transporter"}),                     # pstS
    ("K03407", {"sensor_kinase"}, {"response_regulator"}),                              # cheA
    ("K03194", {"type_iv_secretion"}, {"type_vi_secretion"}),                           # virB1
    ("K17836", {"beta_lactamase"}, {"aminoglycoside_resistance"}),                      # penP
])
def test_corrected_entries(ko, present, absent):
    preds = set(KEGG_TO_PREDICATES[ko])
    assert present <= preds
    assert not absent & preds


@pytest.mark.parametrize("ko", ["K15826", "K14129", "K14130"])
def test_non_hydrogenase_kos_are_unmapped(ko):
    assert ko not in KEGG_TO_PREDICATES


@pytest.mark.parametrize(("ec", "present", "absent"), [
    ("2.1.2.1", {"transferase"}, {"methyltransferase"}),
    ("2.1.1.37", {"methyltransferase"}, set()),
    ("3.6.5.2", {"gtp_binding"}, {"atpase"}),
    ("3.6.1.1", {"hydrolase"}, {"atpase", "atp_binding"}),
    ("4.1.3.1", {"lyase"}, {"decarboxylase"}),
    ("7.6.2.1", {"atp_binding"}, {"gtp_binding"}),
    ("2.7.8.5", {"transferase"}, {"kinase"}),
    ("1.4.1.3", {"oxidoreductase"}, {"aminotransferase"}),
    ("2.6.1.1", {"aminotransferase", "plp_binding"}, set()),
])
def test_ec_classes(ec, present, absent):
    preds = set(get_predicates_for_ec(ec))
    assert present <= preds
    assert not absent & preds
