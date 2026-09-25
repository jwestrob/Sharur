# Maintainer shortcuts. SNAPSHOT is a directory from `make snapshots`.
SNAPSHOT ?= data/reference/snapshots/$(shell date +%Y-%m-%d)
PYTHON ?= python

.PHONY: test lint snapshots snapshots-with-links predicate-maps recount kegg

test:
	$(PYTHON) -m pytest --override-ini addopts="" -q

lint:
	$(PYTHON) -m ruff check --select E9,F63,F7,F82 sharur src/ingest tests

# Reference inputs for the predicate maps (Pfam, GO, ENZYME, Swiss-Prot, KEGG REST).
snapshots:
	$(PYTHON) scripts/fetch_reference_snapshots.py --dest $(SNAPSHOT)

# Also KEGG gene -> KO links for Swiss-Prot organisms (~2,200 KEGG REST requests).
snapshots-with-links:
	$(PYTHON) scripts/fetch_reference_snapshots.py --dest $(SNAPSHOT) --kegg-links

# Rebuild the HydDB x KOfam snapshot, KO consensus, local KEGG map, Pfam and CAZy maps.
predicate-maps:
	$(PYTHON) scripts/rebuild_predicate_maps.py --snapshot $(SNAPSHOT)

# Re-verify every map pair and recount reviewed-protein consensus from the snapshot.
recount:
	SHARUR_SWISSPROT=$(SNAPSHOT)/uniprot_sprot.dat.gz SHARUR_GO_OBO=$(SNAPSHOT)/go-basic.obo \
	SHARUR_SWISSPROT_KEGG=$(SNAPSHOT)/swissprot_kegg_ko.tsv \
	$(PYTHON) -m pytest --override-ini addopts="" -q tests/test_pfam_map_integrity.py \
		tests/test_kegg_map_integrity.py tests/test_cazy_map_integrity.py tests/test_kegg_rules.py

# User setup: build the KEGG map on this machine (academic KEGG REST use).
kegg:
	sharur setup-kegg
