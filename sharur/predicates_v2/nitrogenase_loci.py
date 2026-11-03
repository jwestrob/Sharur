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
   ``WINDOW`` genes downstream of it, before the next Fe-protein KO gene on
   that strand (a tandem nifH pair calls one gene set, from the nifH nearer the
   catalytic genes).
3. Alpha slot KO K02586 (nifD) or K22896 (vnfD); beta slot KO K02591 (nifK) or
   K22897 (vnfK). In the methanogen-type operon (nifH-nifI1-nifI2-nifD-nifK)
   KOfam's best hit for the beta subunit is the NifN group (K02592), so K02592
   fills the beta slot when a nifI1/nifI2 gene (K02589/K02590) lies between the
   Fe protein and the alpha slot.
4. Each of the three KOs scores at or above its KOfam threshold, on the score
   that KOfam's ``score_type`` names for the KO: the full-sequence bitscore
   (nifH, nifD, the vnf genes) or the best-domain bitscore (nifK, the NifN group).
5. None of the three genes is partial (cut by a contig end). A gene without a
   Prodigal partial flag counts as partial when it is the first or last gene of
   a linear contig, since only those genes can run off a contig end.

Every other Fe-protein KO gene, and every alpha/beta KO gene outside a called
gene set, becomes a ``locus_type = 'nitrogenase_candidate'`` row whose ``call``
states what was observed (see ``CALLS``). A gene set cut by a contig end gets
``contig_edge`` and never ``confirmed``. Each gene takes one structural role
(fe_protein, alpha, beta, anchor or partner) in at most one gene set.

Strands are read in any encoding Sharur datasets carry (``+``/``-``,
``1``/``-1``). On a circular contig a gene set may cross the origin: its
members are listed in contig order across the origin (..., n-1, 0, 1, ...), its
``locus_id`` names the first and last member index in that order (so the first
index is the larger one), and its ``start``/``end_coord`` cover the members on
the Fe protein's side of the origin (for an unanchored set, its 5'-most gene),
because every ``loci`` reader expects ``start <= end_coord``. Such rows carry
``wraps_origin: true`` and both arcs as ``segments`` in their metadata.

The KO score: a KO loader records which bitscore its rows store
(:mod:`sharur.kofam`). Rows loaded before that record existed store one domain's
bitscore (Stage 07) or the full-sequence bitscore (other loaders), and the score
KOfam names for a ``full`` KO, or for a ``domain`` KO, may be absent from the
database. The threshold check then compares the stored score and flags the gene
set: its slot members carry ``ko_score_exact: false`` and the locus
``slot_ko_scores_exact: false``; the call itself stands. Astra and Aksha write a
KOFAM hit only when its full-sequence score and the domain's score both reach the
profile's embedded threshold (``ko_list``'s value set as GA1 and GA2), so on rows
from those searches the two scores agree about a threshold from the same
``ko_list``. Reloading the KO rows replaces the flag with KOfam's own comparison.

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
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any

from sharur.kofam import KoThreshold, compares_as_kofam, ko_score_bases, read_ko_list
from sharur.predicates.pfam_identity import normalize_pfam_accession
from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet


if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping
    from pathlib import Path


logger = logging.getLogger(__name__)

CALLER = "sharur.nitrogenase_loci"
CALLER_VERSION = "2"  # 2: KOfam score_type picks the full-sequence or best-domain score
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

_STRAND_ALIASES = {"+": "+", "+1": "+", "1": "+", "forward": "+", "f": "+",
                   "-": "-", "-1": "-", "reverse": "-", "r": "-"}


def normalize_strand(value: Any) -> str:
    """``+`` or ``-`` for every strand encoding Sharur datasets carry; ``?`` for anything else.

    Accepts the values :func:`sharur.operators.cases._normalize_strand` accepts. That
    helper cannot be imported here: ``sharur.operators`` imports the V2 persistence
    layer, which imports this module. A ``?`` gene matches only other ``?`` genes.
    """
    if value is None:
        return "?"
    return _STRAND_ALIASES.get(str(value).strip().lower(), "?")


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
    partial_inferred: bool = False  # no Prodigal flag: ``partial`` comes from the gene's contig position
    ko_score_exact: bool = True  # ``ko_passes`` compared the score KOfam's score_type names for ``ko``

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
    """One called gene set, ready for ``loci`` / ``locus_proteins``.

    ``members`` run in contig order; on a circular contig that order crosses the
    origin when the gene set does (``wraps_origin``).
    """

    contig_id: str
    call: str
    strand: str
    members: list[tuple[Gene, str]]  # (gene, role) in contig order
    evidence: dict[str, Any] = field(default_factory=dict)

    @property
    def locus_type(self) -> str:
        return CONFIRMED_TYPE if self.call == "confirmed" else CANDIDATE_TYPE

    @property
    def wraps_origin(self) -> bool:
        indexes = [gene.gene_index for gene, _ in self.members]
        return any(b < a for a, b in zip(indexes, indexes[1:]))

    @property
    def locus_id(self) -> str:
        first, last = self.members[0][0].gene_index, self.members[-1][0].gene_index
        return f"{self.locus_type}:{self.contig_id}:{first}-{last}"

    def segments(self) -> list[tuple[int, int]]:
        """(start, end) of each stretch of members between contig ends, in member order."""
        runs: list[list[Gene]] = [[]]
        for gene, _ in self.members:
            if runs[-1] and gene.gene_index < runs[-1][-1].gene_index:
                runs.append([])
            runs[-1].append(gene)
        return [(min(g.start for g in run), max(g.end for g in run)) for run in runs]

    def span(self) -> tuple[int, int]:
        """``start``/``end_coord`` for the loci row (always start <= end).

        A gene set that crosses a circular contig's origin is stored as the stretch
        holding its 5'-most member (the Fe protein of an anchored set).
        """
        segments = self.segments()
        if len(segments) == 1:
            return segments[0]
        head = self.members[-1][0] if self.strand == "-" else self.members[0][0]
        return next(seg for seg in segments if seg[0] <= head.start and head.end <= seg[1])

    def metadata(self) -> dict[str, Any]:
        members = []
        for gene, role in self.members:
            member = {"protein_id": gene.protein_id, "role": role, "ko": gene.ko,
                      "ko_passes": gene.ko_passes, "length": gene.length, "partial": gene.partial}
            if gene.partial_inferred and gene.partial:
                member["partial_inferred"] = True  # no Prodigal flag; first or last gene of a linear contig
            if gene.ko is not None and not gene.ko_score_exact:
                member["ko_score_exact"] = False  # the stored score may differ from the one KOfam compares
            members.append(member)
        metadata = {
            "caller": CALLER,
            "caller_version": CALLER_VERSION,
            "call": self.call,
            "call_reason": CALLS[self.call],
            "window": WINDOW,
            "strand": self.strand,
            "members": members,
            **self.evidence,
        }
        if self.wraps_origin:
            metadata["wraps_origin"] = True
            metadata["segments"] = [list(seg) for seg in self.segments()]
        return metadata


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

    def walk(self, origin: Gene, step: int, count: int, stop=None) -> tuple[list[Gene], bool]:
        """Up to ``count`` genes from ``origin`` in direction ``step``; True if a contig end cuts the walk.

        The walk ends before the first gene for which ``stop(gene)`` holds. On a
        circular contig it visits each other gene at most once.
        """
        if self.contig.circular:
            count = min(count, max(self.contig.n_genes - 1, 0))
        genes, cut = [], False
        for offset in range(1, count + 1):
            gene, past_end = self.at(origin.gene_index + step * offset)
            if past_end:
                cut = True
                break
            if gene is None or gene.protein_id == origin.protein_id:
                continue
            if stop is not None and stop(gene):
                break
            genes.append(gene)
        return genes, cut

    def stretch(self, first: int, last: int) -> list[Gene]:
        """Genes from index ``first`` to ``last`` in contig order, across the origin of a circular contig."""
        n = self.contig.n_genes
        if self.contig.circular and n > 0:
            indexes = [(first + k) % n for k in range((last - first) % n + 1)]
        else:
            indexes = list(range(first, last + 1))
        return [self.by_index[i] for i in indexes if i in self.by_index]

    def distance(self, lower: Gene, upper: Gene) -> int:
        """Genes from ``lower`` forward to ``upper`` (across the origin of a circular contig)."""
        gap = upper.gene_index - lower.gene_index
        return gap % self.contig.n_genes if self.contig.circular and self.contig.n_genes > 0 else gap


def _direction(gene: Gene) -> int:
    return -1 if gene.strand == "-" else 1


def _call_anchor(contig: _ContigGenes, anchor: Gene) -> Locus:
    step = _direction(anchor)
    # The next Fe-protein gene on this strand starts its own gene set.
    downstream, cut = contig.walk(anchor, step, WINDOW,
                                  stop=lambda gene: gene.strand == anchor.strand and gene.ko in FE_PROTEIN_KOS)
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

    tail = beta or alpha or anchor
    first, last = (anchor, tail) if step == 1 else (tail, anchor)
    confirmed = call == "confirmed"
    roles = {anchor.protein_id: "fe_protein" if confirmed else "anchor"}
    if alpha is not None:
        roles[alpha.protein_id] = "alpha" if confirmed else "partner"
    if beta is not None:
        roles[beta.protein_id] = "beta" if confirmed else "partner"
    members = []
    for gene in contig.stretch(first.gene_index, last.gene_index):
        role = roles.get(gene.protein_id)
        if role is None:
            role = "nifI" if gene.ko in NIFI_KOS else "delta" if confirmed and gene.ko in DELTA_KOS else "other"
        members.append((gene, role))
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
        "slot_ko_scores_exact": all(gene.ko_score_exact for gene in slots if gene.ko is not None),
    }
    return Locus(anchor.contig_id, call, anchor.strand, members, evidence)


def _call_unanchored(contig: _ContigGenes, genes: list[Gene], claimed: set[str] = frozenset()) -> Locus:
    """Alpha/beta KO genes that no Fe-protein gene claimed: one gene set per strand run.

    ``genes`` are in contig order (across the origin of a circular contig). A run
    holding a single alpha/beta gene takes the next Oxidored_nitro-family gene
    downstream on its strand as its partner (in methanogen-type operons the beta
    subunit carries the NifN-group KO, which seeds nothing on its own), unless a
    called gene set already holds that gene (``claimed``).
    """
    step = _direction(genes[0])
    head, tail = (genes[0], genes[-1]) if step == 1 else (genes[-1], genes[0])  # 5'-most, 3'-most
    run = list(genes)
    if len(run) == 1:
        downstream, _ = contig.walk(tail, step, WINDOW)
        partner = next((gene for gene in downstream if gene.strand == tail.strand and gene.catalytic_family), None)
        if partner is not None and partner.protein_id not in claimed:
            run.append(partner)
            tail = partner
    upstream, cut = contig.walk(head, -step, WINDOW)
    call = "contig_edge" if cut or any(gene.partial for gene in run) else "no_fe_protein"
    first, last = (head, tail) if step == 1 else (tail, head)
    partners = {gene.protein_id for gene in run}
    members = []
    for gene in contig.stretch(first.gene_index, last.gene_index):
        role = "partner" if gene.protein_id in partners else "nifI" if gene.ko in NIFI_KOS else "other"
        members.append((gene, role))
    evidence = {
        "fe_protein_upstream": any(gene.ko in FE_PROTEIN_KOS for gene in upstream),
        "window_cut_by_contig_end": cut,
    }
    return Locus(head.contig_id, call, head.strand, members, evidence)


def _strand_runs(contig: _ContigGenes, genes: list[Gene]) -> list[list[Gene]]:
    """Group genes (in contig order) into runs: same strand, each within ``WINDOW`` genes of the previous.

    Genes on the other strand do not break a run. On a circular contig a run that
    ends near the last gene joins a same-strand run that starts near gene 0.
    """
    runs: list[list[Gene]] = []
    open_runs: dict[str, list[Gene]] = {}
    for gene in genes:
        run = open_runs.get(gene.strand)
        if run is not None and contig.distance(run[-1], gene) <= WINDOW:
            run.append(gene)
        else:
            run = [gene]
            runs.append(run)
            open_runs[gene.strand] = run
    if contig.contig.circular:
        for strand in {run[0].strand for run in runs}:
            same = [run for run in runs if run[0].strand == strand]
            if len(same) > 1 and contig.distance(same[-1][-1], same[0][0]) <= WINDOW:
                same[-1].extend(same[0])
                runs = [run for run in runs if run is not same[0]]
    return runs


def call_loci(contigs: Mapping[str, Contig], genes: Iterable[Gene]) -> list[Locus]:
    """Call every Fe-protein KO gene and every unclaimed alpha/beta KO gene."""
    by_contig: dict[str, list[Gene]] = defaultdict(list)
    for gene in genes:
        strand = normalize_strand(gene.strand)
        by_contig[gene.contig_id].append(gene if strand == gene.strand else replace(gene, strand=strand))

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

        seeds = [gene for gene in contig_genes
                 if (gene.ko in ALPHA_KOS or gene.ko in BETA_KOS) and gene.protein_id not in claimed]
        loci.extend(_call_unanchored(contig, run, claimed) for run in _strand_runs(contig, seeds))

    # One gene set per id (a safeguard: anchored gene sets end before the next same-strand Fe protein).
    unique: dict[str, Locus] = {}
    for locus in loci:
        unique.setdefault(locus.locus_id, locus)
    return list(unique.values())


# --------------------------------------------------------------------------- #
# Database I/O
# --------------------------------------------------------------------------- #


def _kofam_thresholds(kos: Iterable[str]) -> dict[str, KoThreshold] | None:
    """KOfam thresholds and score types for ``kos``; None when ko_list is absent."""
    from sharur.health import kofam_ko_list  # noqa: PLC0415

    path: Path | None = kofam_ko_list()
    return None if path is None else read_ko_list(path, kos)


def _ko_check(ko: str, score: float | None, evalue: float | None,
              thresholds: dict[str, KoThreshold] | None, basis: str | None) -> tuple[bool, bool]:
    """Whether a stored KO call passes KOfam's threshold, and whether that compared KOfam's own score.

    ``basis`` is the score the KO row's loader recorded storing (:mod:`sharur.kofam`);
    the comparison is KOfam's when that is the score the KO's ``score_type`` names.
    KOs without a threshold are judged on the full-sequence E-value, which every
    loader stores.
    """
    if thresholds is None:
        return True, False  # no ko_list: the stored KO call stands (recorded in metadata)
    entry = thresholds.get(ko)
    if entry is None or entry.threshold is None:
        return evalue is not None and evalue <= KOFAM_NO_THRESHOLD_EVALUE, True
    return score is not None and score >= entry.threshold, compares_as_kofam(basis, entry.score_type)


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
    bases = ko_score_bases(conn)
    best_ko: dict[str, tuple[bool, float, str, bool]] = {}
    pfams: dict[str, set[str]] = defaultdict(set)
    for protein_id, source, accession, name, score, evalue in hits:
        if source == "pfam":
            for key in (normalize_pfam_accession(accession), (name or "").strip()):
                if key in CONTEXT_PFAMS:
                    pfams[protein_id].add(key)
            continue
        passes, exact = _ko_check(accession, score, evalue, thresholds, bases.get(source))
        candidate = (passes, score if score is not None else float("-inf"), accession, exact)
        if protein_id not in best_ko or candidate > best_ko[protein_id]:
            best_ko[protein_id] = candidate
    genes = []
    for protein_id, contig_id, gene_index, strand, start, end, length, partial in proteins:
        passes, _, ko, exact = best_ko.get(protein_id, (False, 0.0, None, True))
        index = int(gene_index)
        if partial is None or str(partial).strip() == "":
            # No Prodigal flag: only the first or last gene of a linear contig can run off its end.
            contig = contigs[contig_id]
            is_partial = not contig.circular and index in (0, contig.n_genes - 1)
            inferred = True
        else:
            is_partial = bool(partial) and str(partial).strip("0") != ""
            inferred = False
        genes.append(Gene(
            protein_id=protein_id, contig_id=contig_id, gene_index=index,
            strand=normalize_strand(strand), start=int(start), end=int(end),
            length=int(length) if length is not None else None,
            partial=is_partial, partial_inferred=inferred,
            ko=ko, ko_passes=passes, pfams=frozenset(pfams.get(protein_id, ())), ko_score_exact=exact,
        ))
    return contigs, genes, thresholds is not None


def _member_states(rows: Iterable[tuple[str, str, Any]]) -> dict[str, set[tuple[str, str, str, str]]]:
    """protein_id -> {(locus_id, locus_type, call, role)}: everything :func:`locus_atoms` reads."""
    states: dict[str, set[tuple[str, str, str, str]]] = defaultdict(set)
    for locus_id, locus_type, metadata in rows:
        if isinstance(metadata, str):
            metadata = json.loads(metadata)
        metadata = metadata or {}
        for member in metadata.get("members", ()):
            states[member.get("protein_id")].add((locus_id, locus_type, metadata.get("call"), member.get("role")))
    return states


def write_loci(conn, loci: list[Locus], *, thresholds_checked: bool = True) -> dict[str, Any]:
    """Replace every nitrogenase locus row; return counts by locus type and call.

    ``summary["changed_proteins"]`` lists the proteins whose locus atoms change with
    the rewrite (a locus, call or role gained or lost); a subset generation must
    regenerate them along with the proteins it was asked for.
    """
    loci_rows, member_rows = [], []
    for locus in sorted(loci, key=lambda item: item.locus_id):
        metadata = locus.metadata()
        metadata["kofam_thresholds_checked"] = thresholds_checked
        start, end = locus.span()
        loci_rows.append((
            locus.locus_id, locus.locus_type, locus.contig_id, start, end,
            1.0 if locus.call == "confirmed" else None, json.dumps(metadata, sort_keys=True),
        ))
        member_rows.extend((locus.locus_id, gene.protein_id, position)
                           for position, (gene, _) in enumerate(locus.members))

    types = ", ".join(f"'{t}'" for t in LOCUS_TYPES)
    by_call: dict[str, int] = defaultdict(int)
    for locus in loci:
        by_call[locus.call] += 1
    # Threshold-gated calls whose slot scores may differ from the ones KOfam compares.
    inexact = sum(1 for locus in loci if locus.call in ("confirmed", "ko_below_threshold")
                  and locus.evidence.get("slot_ko_scores_exact") is False)
    summary = {"loci": len(loci_rows), "members": len(member_rows), "by_call": dict(sorted(by_call.items())),
               "kofam_thresholds_checked": thresholds_checked, "ko_score_inexact_loci": inexact,
               "written": False, "changed_proteins": []}
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

    before = _member_states((row[0], row[1], row[6]) for row in stored_loci)
    after = _member_states((row[0], row[1], row[6]) for row in loci_rows)
    summary["changed_proteins"] = sorted(
        pid for pid in before.keys() | after.keys() if pid is not None and before.get(pid) != after.get(pid))

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
    whether they were, and ``summary["changed_proteins"]`` which proteins' locus
    atoms the rewrite changed.
    """
    conn = store.conn
    tables = {row[0] for row in conn.execute("SHOW TABLES").fetchall()}
    if not {"proteins", "annotations", "loci", "locus_proteins"} <= tables:
        return {"loci": 0, "members": 0, "by_call": {}, "kofam_thresholds_checked": False,
                "ko_score_inexact_loci": 0, "written": False, "changed_proteins": []}
    contigs, genes, checked = fetch_genes(conn)
    if not checked:
        logger.warning("KOfam ko_list not found: nitrogenase slot KOs are taken as called, unthresholded")
    loci = call_loci(contigs, genes)
    summary = write_loci(conn, loci, thresholds_checked=checked)
    if checked and summary["ko_score_inexact_loci"]:
        logger.warning(
            "%d nitrogenase gene sets were thresholded on a stored KO score whose basis the KO loader did not "
            "record (KOfam compares the full-sequence or best-domain bitscore per score_type); reload the KO "
            "rows to compare KOfam's own score. Flagged as slot_ko_scores_exact: false.",
            summary["ko_score_inexact_loci"])
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
    "normalize_strand",
    "refresh_nitrogenase_loci",
    "write_loci",
]
