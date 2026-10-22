"""Nitrogenase structural-locus calls from gene order.

A nitrogenase is a multi-gene system: the Fe protein (NifH, VnfH) and the alpha
(NifD, VnfD) and beta (NifK, VnfK) subunits of the catalytic component. The
NifH and NifD/NifK families also hold paralogs that make no nitrogenase: CfbC/CfbD
(coenzyme F430 synthesis), the Bch/Chl chlorophyllide reductases, and NifH-like
proteins paired with shorter or KO-unassigned NifD-like partners. A single gene's
ortholog call therefore cannot name the enzyme, so this caller reads gene order
and writes one ``loci`` row per gene set:

``locus_type = 'nitrogenase'`` (call ``confirmed``) when all of these hold:

1. Fe protein: a gene whose KO is K02588 (nifH) or K22899 (vnfH).
2. Catalytic slots: the first two Oxidored_nitro-family genes (Pfam
   Oxidored_nitro, or a KO of that family) on the Fe protein's strand within
   ``WINDOW`` genes downstream of it.
3. Alpha slot KO K02586 (nifD) or K22896 (vnfD); beta slot KO K02591 (nifK) or
   K22897 (vnfK). In the methanogen-type operon (nifH-nifI1-nifI2-nifD-nifK)
   KOfam's best hit for the beta subunit is the NifN group (K02592), so K02592
   fills the beta slot when a nifI1/nifI2 gene (K02589/K02590) lies between the
   Fe protein and the alpha slot.
4. Each of the three KOs scores at or above its KOfam threshold.
5. None of the three genes is partial (cut by a contig end).

Every other Fe-protein KO gene, and every alpha/beta KO gene outside a called
gene set, becomes a ``locus_type = 'nitrogenase_candidate'`` row whose ``call``
states what was observed (see ``CALLS``). A gene set cut by a contig end gets
``contig_edge`` and never ``confirmed``.

The window: in the SR-VP archaeal and bacterial snapshots every KO-assigned
alpha and beta partner of a nifH lies 1-5 genes downstream on its strand (alpha
+1 to +3, beta +2 to +5; the nifI1/nifI2 pair accounts for +2); none lies
upstream or on the other strand. ``WINDOW = 6`` keeps one gene of slack, as the
hydrogenase neighborhood check does.

V2 predicate generation reads the rows back (:func:`locus_atoms`): structural
members of a confirmed gene set carry the system-level atom
``nitrogenase_locus``; the anchor and partner genes of a ``contig_edge`` set
carry the quality flag ``nitrogenase_contig_edge``.
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from sharur.predicates.pfam_identity import normalize_pfam_accession
from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet


if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping
    from pathlib import Path


logger = logging.getLogger(__name__)

CALLER = "sharur.nitrogenase_loci"
CALLER_VERSION = "1"
WINDOW = 6
SHORT_PARTNER_AA = 430  # labels uncalled pairs; KO-assigned alpha/beta partners here are >= 439 aa
SOURCE_DB = "nitrogenase_locus"
CONFIRMED_TYPE = "nitrogenase"
CANDIDATE_TYPE = "nitrogenase_candidate"
LOCUS_TYPES = (CONFIRMED_TYPE, CANDIDATE_TYPE)
SYSTEM_ATOM = "nitrogenase_locus"
EDGE_ATOM = "nitrogenase_contig_edge"

FE_PROTEIN_KOS = {"K02588": "nifH", "K22899": "vnfH"}
ALPHA_KOS = {"K02586": "nifD", "K22896": "vnfD"}
BETA_KOS = {"K02591": "nifK", "K22897": "vnfK"}
GROUP_II_BETA_KO = "K02592"  # nifN ortholog group; beta slot only beside nifI1/nifI2
DELTA_KOS = {"K00531": "anfG", "K22898": "vnfG"}
NIFI_KOS = {"K02589": "nifI1", "K02590": "nifI2"}
COFACTOR_KOS = {"K02587": "nifE", "K02592": "nifN", "K22903": "vnfE"}
NIFB_KO = "K02585"
PARALOG_KOS = {
    "K21610": "cfbC", "K21611": "cfbD",
    "K04037": "chlL", "K04038": "chlN", "K04039": "chlB",
    "K11333": "bchX", "K11334": "bchY", "K11335": "bchZ",
}
CATALYTIC_FAMILY_KOS = frozenset(
    {*ALPHA_KOS, *BETA_KOS, *COFACTOR_KOS, "K21611", "K04038", "K04039", "K11334", "K11335"}
)
SEED_KOS = frozenset({*FE_PROTEIN_KOS, *ALPHA_KOS, *BETA_KOS})
CONTEXT_KOS = frozenset(
    {*SEED_KOS, GROUP_II_BETA_KO, *DELTA_KOS, *NIFI_KOS, *COFACTOR_KOS, NIFB_KO, *PARALOG_KOS}
)
CATALYTIC_PFAM = frozenset({"PF00148", "Oxidored_nitro"})
NIFX_PFAM = frozenset({"PF02579", "Nitro_FeMo-Co"})
CONTEXT_PFAMS = CATALYTIC_PFAM | NIFX_PFAM

STRUCTURAL_ROLES = frozenset({"fe_protein", "alpha", "beta", "delta"})
EDGE_FLAG_ROLES = frozenset({"anchor", "partner"})
KOFAM_NO_THRESHOLD_EVALUE = 1e-15  # KOfam's cascade for KOs without a score threshold

CALLS = {
    "confirmed": "Fe protein, alpha and beta slot KOs in operon order, all at KOfam threshold",
    "contig_edge": "gene set cut by a contig end (window runs off the contig or a slot gene is partial)",
    "ko_below_threshold": "nitrogenase gene order, with a slot KO below its KOfam threshold",
    "beta_slot_mismatch": "Fe protein and alpha slot KOs; the beta slot holds another gene or is empty",
    "paralog_ko": "a catalytic slot carries a CfbD or Bch/Chl reductase KO",
    "cofactor_ko_partners": "catalytic slots carry NifE/NifN/VnfE KOs instead of alpha/beta",
    "short_partners": f"two KO-unassigned Oxidored_nitro partners shorter than {SHORT_PARTNER_AA} aa",
    "unassigned_partners": f"two KO-unassigned Oxidored_nitro partners of {SHORT_PARTNER_AA} aa or more",
    "single_partner": "one Oxidored_nitro partner downstream",
    "upstream_partners": "no partner downstream; two or more Oxidored_nitro genes upstream on the same strand",
    "no_partner": "no Oxidored_nitro partner within the window",
    "no_fe_protein": "alpha/beta KO genes with no Fe-protein KO gene upstream within the window",
}


@dataclass(frozen=True)
class Gene:
    """One gene and the annotations the caller reads."""

    protein_id: str
    contig_id: str
    gene_index: int
    strand: str
    start: int
    end: int
    length: int | None = None
    partial: bool = False
    ko: str | None = None
    ko_passes: bool = False
    pfams: frozenset[str] = frozenset()

    @property
    def catalytic_family(self) -> bool:
        return self.ko in CATALYTIC_FAMILY_KOS or bool(self.pfams & CATALYTIC_PFAM)


@dataclass(frozen=True)
class Contig:
    contig_id: str
    n_genes: int  # max gene_index + 1
    circular: bool = False


@dataclass
class Locus:
    """One called gene set, ready for ``loci`` / ``locus_proteins``."""

    contig_id: str
    call: str
    strand: str
    members: list[tuple[Gene, str]]  # (gene, role) in contig order
    evidence: dict[str, Any] = field(default_factory=dict)

    @property
    def locus_type(self) -> str:
        return CONFIRMED_TYPE if self.call == "confirmed" else CANDIDATE_TYPE

    @property
    def locus_id(self) -> str:
        indexes = [gene.gene_index for gene, _ in self.members]
        return f"{self.locus_type}:{self.contig_id}:{min(indexes)}-{max(indexes)}"

    def metadata(self) -> dict[str, Any]:
        return {
            "caller": CALLER,
            "caller_version": CALLER_VERSION,
            "call": self.call,
            "call_reason": CALLS[self.call],
            "window": WINDOW,
            "strand": self.strand,
            "members": [
                {"protein_id": gene.protein_id, "role": role, "ko": gene.ko,
                 "ko_passes": gene.ko_passes, "length": gene.length, "partial": gene.partial}
                for gene, role in self.members
            ],
            **self.evidence,
        }


# --------------------------------------------------------------------------- #
# Calling
# --------------------------------------------------------------------------- #


class _ContigGenes:
    """Gene lookup on one contig, wrapping around circular contigs."""

    def __init__(self, contig: Contig, genes: Iterable[Gene]):
        self.contig = contig
        self.by_index = {gene.gene_index: gene for gene in genes}

    def at(self, index: int) -> tuple[Gene | None, bool]:
        """Gene at ``index`` and whether the index lies past a contig end."""
        n = self.contig.n_genes
        if self.contig.circular and n > 0:
            index %= n
        elif index < 0 or index >= n:
            return None, True
        return self.by_index.get(index), False

    def walk(self, origin: Gene, step: int, count: int) -> tuple[list[Gene], bool]:
        """Up to ``count`` genes from ``origin`` in direction ``step``; True if a contig end cuts the walk."""
        genes, cut = [], False
        for offset in range(1, count + 1):
            gene, past_end = self.at(origin.gene_index + step * offset)
            if past_end:
                cut = True
                break
            if gene is not None and gene.protein_id != origin.protein_id:
                genes.append(gene)
        return genes, cut


def _direction(gene: Gene) -> int:
    return -1 if gene.strand == "-" else 1


def _call_anchor(contig: _ContigGenes, anchor: Gene) -> Locus:
    step = _direction(anchor)
    downstream, cut = contig.walk(anchor, step, WINDOW)
    same_strand = [gene for gene in downstream if gene.strand == anchor.strand]
    partners = [gene for gene in same_strand if gene.catalytic_family]
    alpha = partners[0] if partners else None
    beta = partners[1] if len(partners) > 1 else None
    before_alpha = same_strand[: same_strand.index(alpha)] if alpha else same_strand
    nifi = [gene for gene in before_alpha if gene.ko in NIFI_KOS]

    upstream, _ = contig.walk(anchor, -step, WINDOW)
    upstream_partners = [gene for gene in upstream if gene.strand == anchor.strand and gene.catalytic_family]
    context, _ = contig.walk(anchor, step, WINDOW + 4)
    nifb_nifx = any(gene.ko == NIFB_KO or gene.pfams & NIFX_PFAM for gene in [*upstream, *context])

    beta_ko_fits = beta is not None and (
        beta.ko in BETA_KOS or (beta.ko == GROUP_II_BETA_KO and bool(nifi))
    )
    slots = [gene for gene in (anchor, alpha, beta) if gene is not None]
    partial = any(gene.partial for gene in slots)

    if alpha is not None and alpha.ko in ALPHA_KOS and beta_ko_fits:
        if partial:
            call = "contig_edge"
        elif all(gene.ko_passes for gene in slots):
            call = "confirmed"
        else:
            call = "ko_below_threshold"
    elif partial or (len(partners) < 2 and cut):
        call = "contig_edge"
    elif any(gene.ko in PARALOG_KOS for gene in partners[:2]):
        call = "paralog_ko"
    elif alpha is not None and alpha.ko in ALPHA_KOS:
        call = "beta_slot_mismatch"
    elif len(partners) >= 2:
        if any(gene.ko in COFACTOR_KOS for gene in partners[:2]):
            call = "cofactor_ko_partners"
        elif max(gene.length or 0 for gene in partners[:2]) < SHORT_PARTNER_AA:
            call = "short_partners"
        else:
            call = "unassigned_partners"
    elif len(partners) == 1:
        call = "single_partner"
    elif len(upstream_partners) >= 2:
        call = "upstream_partners"
    else:
        call = "no_partner"

    last = beta or alpha
    span = [anchor]
    if last is not None:
        span += downstream[: downstream.index(last) + 1]
    confirmed = call == "confirmed"
    roles = {anchor.protein_id: "fe_protein" if confirmed else "anchor"}
    if alpha is not None:
        roles[alpha.protein_id] = "alpha" if confirmed else "partner"
    if beta is not None:
        roles[beta.protein_id] = "beta" if confirmed else "partner"
    members = []
    for gene in span:
        role = roles.get(gene.protein_id)
        if role is None:
            role = "nifI" if gene.ko in NIFI_KOS else "delta" if confirmed and gene.ko in DELTA_KOS else "other"
        members.append((gene, role))
    members.sort(key=lambda item: item[0].gene_index)
    evidence = {
        "partner_lengths": [gene.length for gene in partners[:2]],
        "partner_kos": [gene.ko for gene in partners[:2]],
        "partners_in_window": len(partners),
        "upstream_partners": len(upstream_partners),
        "upstream_partner_lengths": [gene.length for gene in upstream_partners[:2]],
        "nifI": bool(nifi),
        "nifB_or_nifX": nifb_nifx,
        "window_cut_by_contig_end": cut,
        "group_ii_beta": bool(beta is not None and beta.ko == GROUP_II_BETA_KO and beta_ko_fits),
    }
    return Locus(anchor.contig_id, call, anchor.strand, members, evidence)


def _call_unanchored(contig: _ContigGenes, genes: list[Gene]) -> Locus:
    """Alpha/beta KO genes that no Fe-protein gene claimed: one gene set per strand run.

    A run holding a single alpha/beta gene takes the next Oxidored_nitro-family gene
    downstream on its strand as its partner (in methanogen-type operons the beta
    subunit carries the NifN-group KO, which seeds nothing on its own).
    """
    step = _direction(genes[0])
    first, last = (genes[0], genes[-1]) if step == 1 else (genes[-1], genes[0])
    run = list(genes)
    if len(run) == 1:
        downstream, _ = contig.walk(last, step, WINDOW)
        partner = next((gene for gene in downstream if gene.strand == last.strand and gene.catalytic_family), None)
        if partner is not None:
            run.append(partner)
    upstream, cut = contig.walk(first, -step, WINDOW)
    call = "contig_edge" if cut or any(gene.partial for gene in run) else "no_fe_protein"
    lo = min(gene.gene_index for gene in run)
    hi = max(gene.gene_index for gene in run)
    partners = {gene.protein_id for gene in run}
    members = []
    for index in range(lo, hi + 1):
        gene, _ = contig.at(index)
        if gene is None:
            continue
        role = "partner" if gene.protein_id in partners else "nifI" if gene.ko in NIFI_KOS else "other"
        members.append((gene, role))
    evidence = {
        "fe_protein_upstream": any(gene.ko in FE_PROTEIN_KOS for gene in upstream),
        "window_cut_by_contig_end": cut,
    }
    return Locus(first.contig_id, call, first.strand, members, evidence)


def call_loci(contigs: Mapping[str, Contig], genes: Iterable[Gene]) -> list[Locus]:
    """Call every Fe-protein KO gene and every unclaimed alpha/beta KO gene."""
    by_contig: dict[str, list[Gene]] = defaultdict(list)
    for gene in genes:
        by_contig[gene.contig_id].append(gene)

    loci: list[Locus] = []
    for contig_id in sorted(by_contig):
        contig_genes = sorted(by_contig[contig_id], key=lambda gene: gene.gene_index)
        contig = _ContigGenes(contigs[contig_id], contig_genes)
        claimed: set[str] = set()
        for anchor in contig_genes:
            if anchor.ko not in FE_PROTEIN_KOS:
                continue
            locus = _call_anchor(contig, anchor)
            loci.append(locus)
            claimed.update(gene.protein_id for gene, role in locus.members if role in ("alpha", "beta", "partner"))

        runs: list[list[Gene]] = []
        for gene in contig_genes:
            if (gene.ko not in ALPHA_KOS and gene.ko not in BETA_KOS) or gene.protein_id in claimed:
                continue
            last = runs[-1][-1] if runs else None
            if last is not None and last.strand == gene.strand and gene.gene_index - last.gene_index <= WINDOW:
                runs[-1].append(gene)
            else:
                runs.append([gene])
        loci.extend(_call_unanchored(contig, run) for run in runs)

    # One gene set per id: two Fe-protein genes can share a span only when one holds no partner.
    unique: dict[str, Locus] = {}
    for locus in loci:
        unique.setdefault(locus.locus_id, locus)
    return list(unique.values())


# --------------------------------------------------------------------------- #
# Database I/O
# --------------------------------------------------------------------------- #


def _kofam_thresholds(kos: Iterable[str]) -> dict[str, float | None] | None:
    """KOfam score thresholds for ``kos`` (None value: KOfam has none); None when ko_list is absent."""
    from sharur.health import kofam_ko_list  # noqa: PLC0415

    path: Path | None = kofam_ko_list()
    if path is None:
        return None
    wanted = set(kos)
    thresholds: dict[str, float | None] = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        ko_idx, thr_idx = header.index("knum"), header.index("threshold")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[ko_idx] in wanted:
                value = parts[thr_idx].strip()
                thresholds[parts[ko_idx]] = None if value in ("", "-") else float(value)
    return thresholds


def _ko_passes(ko: str, score: float | None, evalue: float | None,
               thresholds: dict[str, float | None] | None) -> bool:
    if thresholds is None:
        return True  # no ko_list: the stored KO call stands (recorded in metadata)
    threshold = thresholds.get(ko)
    if threshold is None:
        return evalue is not None and evalue <= KOFAM_NO_THRESHOLD_EVALUE
    return score is not None and score >= threshold


def fetch_genes(conn) -> tuple[dict[str, Contig], list[Gene], bool]:
    """Genes within ``WINDOW + 4`` of every Nif-family seed, with their KO and Pfam flags.

    Returns contigs, genes, and whether KOfam thresholds were checked.
    """
    seed_kos = sorted(SEED_KOS)
    seeds = conn.execute(
        f"""SELECT DISTINCT p.contig_id, p.gene_index
            FROM annotations a JOIN proteins p USING (protein_id)
            WHERE LOWER(a.source) IN ('kegg', 'kofam')
              AND a.accession IN ({', '.join('?' * len(seed_kos))})
              AND p.gene_index IS NOT NULL""",
        seed_kos,
    ).fetchall()
    if not seeds:
        return {}, [], True
    seed_contigs = sorted({contig_id for contig_id, _ in seeds})
    conn.execute("CREATE OR REPLACE TEMP TABLE _nif_seed_contigs (contig_id VARCHAR)")
    conn.executemany("INSERT INTO _nif_seed_contigs VALUES (?)", [(c,) for c in seed_contigs])
    contigs = {
        contig_id: Contig(contig_id, int(max_index) + 1, bool(circular))
        for contig_id, max_index, circular in conn.execute(
            """SELECT p.contig_id, MAX(p.gene_index), COALESCE(ANY_VALUE(c.is_circular), FALSE)
               FROM proteins p JOIN _nif_seed_contigs s USING (contig_id)
               LEFT JOIN contigs c ON c.contig_id = p.contig_id
               GROUP BY p.contig_id"""
        ).fetchall()
    }
    reach = WINDOW + 4
    wanted: set[tuple[str, int]] = set()
    for contig_id, index in seeds:
        n = contigs[contig_id].n_genes
        for offset in range(-reach, reach + 1):
            j = index + offset
            if contigs[contig_id].circular:
                wanted.add((contig_id, j % n))
            elif 0 <= j < n:
                wanted.add((contig_id, j))
    conn.execute("CREATE OR REPLACE TEMP TABLE _nif_window (contig_id VARCHAR, gene_index INTEGER)")
    conn.executemany("INSERT INTO _nif_window VALUES (?, ?)", sorted(wanted))
    proteins = conn.execute(
        """SELECT p.protein_id, p.contig_id, p.gene_index, p.strand, p.start, p.end_coord,
                  p.sequence_length, p.partial
           FROM proteins p JOIN _nif_window w
             ON w.contig_id = p.contig_id AND w.gene_index = p.gene_index"""
    ).fetchall()
    context_kos = sorted(CONTEXT_KOS)
    context_pfams = sorted(CONTEXT_PFAMS)
    hits = conn.execute(
        f"""SELECT a.protein_id, LOWER(a.source), a.accession, a.name, a.score, a.evalue
            FROM annotations a
            JOIN proteins p USING (protein_id)
            JOIN _nif_window w ON w.contig_id = p.contig_id AND w.gene_index = p.gene_index
            WHERE (LOWER(a.source) IN ('kegg', 'kofam') AND a.accession IN ({', '.join('?' * len(context_kos))}))
               OR (LOWER(a.source) = 'pfam'
                   AND (regexp_replace(a.accession, '\\.[0-9]+$', '') IN ({', '.join('?' * len(context_pfams))})
                        OR a.name IN ({', '.join('?' * len(context_pfams))})))""",
        [*context_kos, *context_pfams, *context_pfams],
    ).fetchall()
    conn.execute("DROP TABLE IF EXISTS _nif_window")
    conn.execute("DROP TABLE IF EXISTS _nif_seed_contigs")

    thresholds = _kofam_thresholds(CONTEXT_KOS)
    best_ko: dict[str, tuple[bool, float, str]] = {}
    pfams: dict[str, set[str]] = defaultdict(set)
    for protein_id, source, accession, name, score, evalue in hits:
        if source == "pfam":
            for key in (normalize_pfam_accession(accession), (name or "").strip()):
                if key in CONTEXT_PFAMS:
                    pfams[protein_id].add(key)
            continue
        passes = _ko_passes(accession, score, evalue, thresholds)
        candidate = (passes, score if score is not None else float("-inf"), accession)
        if protein_id not in best_ko or candidate > best_ko[protein_id]:
            best_ko[protein_id] = candidate
    genes = []
    for protein_id, contig_id, gene_index, strand, start, end, length, partial in proteins:
        passes, _, ko = best_ko.get(protein_id, (False, 0.0, None))
        genes.append(Gene(
            protein_id=protein_id, contig_id=contig_id, gene_index=int(gene_index),
            strand=strand or "+", start=int(start), end=int(end),
            length=int(length) if length is not None else None,
            partial=bool(partial) and str(partial).strip("0") != "",
            ko=ko, ko_passes=passes, pfams=frozenset(pfams.get(protein_id, ())),
        ))
    return contigs, genes, thresholds is not None


def write_loci(conn, loci: list[Locus], *, thresholds_checked: bool = True) -> dict[str, Any]:
    """Replace every nitrogenase locus row; return counts by locus type and call."""
    loci_rows, member_rows = [], []
    for locus in sorted(loci, key=lambda item: item.locus_id):
        genes = [gene for gene, _ in locus.members]
        metadata = locus.metadata()
        metadata["kofam_thresholds_checked"] = thresholds_checked
        loci_rows.append((
            locus.locus_id, locus.locus_type, locus.contig_id,
            min(gene.start for gene in genes), max(gene.end for gene in genes),
            1.0 if locus.call == "confirmed" else None, json.dumps(metadata, sort_keys=True),
        ))
        member_rows.extend((locus.locus_id, gene.protein_id, position)
                           for position, (gene, _) in enumerate(locus.members))

    types = ", ".join(f"'{t}'" for t in LOCUS_TYPES)
    by_call: dict[str, int] = defaultdict(int)
    for locus in loci:
        by_call[locus.call] += 1
    summary = {"loci": len(loci_rows), "members": len(member_rows), "by_call": dict(sorted(by_call.items())),
               "kofam_thresholds_checked": thresholds_checked, "written": False}
    stored_loci = conn.execute(
        f"""SELECT locus_id, locus_type, contig_id, start, end_coord, confidence, CAST(metadata AS VARCHAR)
            FROM loci WHERE locus_type IN ({types}) ORDER BY locus_id"""
    ).fetchall()
    stored_members = conn.execute(
        f"""SELECT lp.locus_id, lp.protein_id, lp.position FROM locus_proteins lp
            JOIN loci l ON l.locus_id = lp.locus_id WHERE l.locus_type IN ({types})"""
    ).fetchall()
    if stored_loci == loci_rows and sorted(stored_members) == sorted(member_rows):
        return summary  # unchanged calls leave the tables (and their storage blocks) untouched

    # DuckDB checks a foreign key against the committed child table, so member rows
    # are deleted and committed before their parent loci (as in sharur.bgc).
    conn.execute(f"DELETE FROM locus_proteins WHERE locus_id IN (SELECT locus_id FROM loci WHERE locus_type IN ({types}))")
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute(f"DELETE FROM loci WHERE locus_type IN ({types})")
        if loci_rows:
            conn.executemany(
                """INSERT INTO loci (locus_id, locus_type, contig_id, start, end_coord, confidence, metadata)
                   VALUES (?, ?, ?, ?, ?, ?, ?)""",
                loci_rows,
            )
            conn.executemany("INSERT INTO locus_proteins (locus_id, protein_id, position) VALUES (?, ?, ?)",
                             member_rows)
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    summary["written"] = True
    return summary


def refresh_nitrogenase_loci(store) -> dict[str, Any]:
    """Call nitrogenase gene sets from the database's annotations and replace their loci rows.

    The rows are rewritten only when the calls change; ``summary["written"]`` says
    whether they were.
    """
    conn = store.conn
    tables = {row[0] for row in conn.execute("SHOW TABLES").fetchall()}
    if not {"proteins", "annotations", "loci", "locus_proteins"} <= tables:
        return {"loci": 0, "members": 0, "by_call": {}, "kofam_thresholds_checked": False, "written": False}
    contigs, genes, checked = fetch_genes(conn)
    if not checked:
        logger.warning("KOfam ko_list not found: nitrogenase slot KOs are taken as called, unthresholded")
    loci = call_loci(contigs, genes)
    summary = write_loci(conn, loci, thresholds_checked=checked)
    logger.info("Nitrogenase loci (%s): %s", "rewritten" if summary["written"] else "unchanged", summary["by_call"])
    return summary


# --------------------------------------------------------------------------- #
# Atoms
# --------------------------------------------------------------------------- #


def locus_atoms(protein_id: str, locus_id: str, locus_type: str, metadata: Any) -> list[SemanticAtom]:
    """Context atoms one locus member carries."""
    if isinstance(metadata, str):
        metadata = json.loads(metadata)
    metadata = metadata or {}
    role = next((m.get("role") for m in metadata.get("members", ()) if m.get("protein_id") == protein_id), None)
    if locus_type == CONFIRMED_TYPE and metadata.get("call") == "confirmed" and role in STRUCTURAL_ROLES:
        return [SemanticAtom(
            protein_id=protein_id, atom_id=SYSTEM_ATOM, facet=SemanticFacet.activity,
            relation=ClaimRelation.implies, source_accession=locus_id, source_db=SOURCE_DB,
        )]
    if locus_type == CANDIDATE_TYPE and metadata.get("call") == "contig_edge" and role in EDGE_FLAG_ROLES:
        return [SemanticAtom(
            protein_id=protein_id, atom_id=EDGE_ATOM, facet=SemanticFacet.quality_flag,
            relation=ClaimRelation.flags, source_accession=locus_id, source_db=SOURCE_DB,
        )]
    return []


__all__ = [
    "CALLS",
    "CANDIDATE_TYPE",
    "CONFIRMED_TYPE",
    "EDGE_ATOM",
    "LOCUS_TYPES",
    "SOURCE_DB",
    "SYSTEM_ATOM",
    "WINDOW",
    "Contig",
    "Gene",
    "Locus",
    "call_loci",
    "fetch_genes",
    "locus_atoms",
    "refresh_nitrogenase_loci",
    "write_loci",
]
