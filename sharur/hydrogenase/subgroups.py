"""HydDB subgroup interpretations and the Sharur predicates each one supports.

Source: Søndergaard D, Pedersen CNS, Greening C (2016) HydDB: a web tool for
hydrogenase classification and analysis. Sci Rep 6:34212, Table 1 and Methods
(doi:10.1038/srep34212; CC BY 4.0). ``reference_role`` quotes Table 1.

The installed reference (HydDB MM2022 FASTA/DIAMOND) also carries labels that
postdate Table 1, such as ``[NiFe] Group_1l``. Those keep their exact label and
receive structural predicates only, with ``status="unverified"``, until their
interpretation is verified against a versioned source.

Predicates are split by kind:

- ``structural``: group membership (``nife_group1``...); the type follows from the hierarchy.
- ``functional``: physiological traits stated for a characterized subgroup.
  Putative and unresolved subgroups emit none. FeFe Group A subtypes emit none
  because HydDB distinguishes A1-A4 by the downstream gene (GltA/GltD -> A2,
  NuoF -> A3, HycB -> A4), which a nearest-reference sequence match lacks.

Every Sharur assignment is a nearest-reference DIAMOND match against HydDB
references, so downstream code treats these predicates as provisional support.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


REFERENCE_SOURCE = "Søndergaard et al. 2016 Sci Rep 6:34212, Table 1"

CHARACTERIZED = "characterized"
PUTATIVE = "putative"
UNRESOLVED = "unresolved"
UNVERIFIED = "unverified"


@dataclass(frozen=True)
class Subgroup:
    hyd_type: str  # NiFe | FeFe | Fe
    subgroup: str  # exact reference label, e.g. Group_1h
    reference_name: str
    reference_role: str
    status: str
    functional: tuple[str, ...] = ()

    @property
    def label(self) -> str:
        return f"[{self.hyd_type}] {self.subgroup}"

    @property
    def structural(self) -> tuple[str, ...]:
        """Group membership; the hydrogenase type follows from the hierarchy."""
        group = self.subgroup.removeprefix("Group_")[:1]
        if self.hyd_type == "NiFe":
            return (f"nife_group{group}",)
        if self.hyd_type == "FeFe":
            return (f"fefe_group{group}",)
        return ("fe_only_hydrogenase",)

    @property
    def predicates(self) -> tuple[str, ...]:
        return self.structural + self.functional


def _s(hyd_type, subgroup, name, role, status, *functional):
    return Subgroup(hyd_type, subgroup, name, role, status, tuple(functional))


_UPTAKE = "uptake_hydrogenase"
_EVOLVING = "h2_evolving"
_ION = "energy_conserving_hydrogenase"
_FD = "ferredoxin_coupled"
_BIDIR = "bidirectional_hydrogenase"

SUBGROUPS: dict[tuple[str, str], Subgroup] = {
    (s.hyd_type, s.subgroup): s
    for s in [
        # [NiFe] Group 1: respiratory H2-uptake
        _s("NiFe", "Group_1a", "Periplasmic",
           "Electron input for sulfate, metal, and organohalide respiration. [NiFeSe] variants.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1b", "Prototypical",
           "Electron input for sulfate, fumarate, metal, and nitrate respiration.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1c", "Hyb-type",
           "Electron input for fumarate, nitrate, and sulfate respiration. Physiologically reversible.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1d", "Oxygen-tolerant",
           "Electron input for aerobic respiration and oxygen-tolerant anaerobic respiration.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1e", "Isp-type",
           "Electron input primarily for sulfur respiration. Physiologically reversible.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1f", "Oxygen-protecting",
           "Unresolved role. May liberate electrons to reduce reactive oxygen species.",
           UNRESOLVED),
        _s("NiFe", "Group_1g", "Crenarchaeota-type",
           "Electron input primarily for sulfur respiration.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1h", "Actinobacteria-type",
           "Electron input for aerobic respiration. Scavenges electrons from atmospheric H2.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1i", "Coriobacteria-type (putative)",
           "Undetermined role. May liberate electrons for anaerobic respiration.",
           PUTATIVE),
        _s("NiFe", "Group_1j", "Archaeoglobi-type",
           "Electron input for sulfate respiration.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1k", "Methanophenazine-reducing",
           "Electron input for methanogenic heterodisulfide respiration.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_1l", "(not in Table 1)",
           "Label present in the installed HydDB MM2022 reference; interpretation unverified.",
           UNVERIFIED),
        # [NiFe] Group 2: alternative and sensory uptake
        _s("NiFe", "Group_2a", "Cyanobacteria-type",
           "Electron input for aerobic respiration. Recycles H2 produced by other cellular processes.",
           CHARACTERIZED, _UPTAKE),
        _s("NiFe", "Group_2b", "Histidine kinase-linked",
           "H2 sensing. Activates two-component system controlling hydrogenase expression.",
           CHARACTERIZED, "h2_sensor"),
        _s("NiFe", "Group_2c", "Diguanylate cyclase-linked (putative)",
           "Undetermined role. May sense H2 and regulate processes through cyclic di-GMP production.",
           PUTATIVE),
        _s("NiFe", "Group_2d", "Aquificae-type",
           "Unresolved role. May generate reductant for carbon fixation or have a regulatory role.",
           UNRESOLVED),
        _s("NiFe", "Group_2e", "Metallosphaera-type (putative)",
           "Undetermined role. May liberate electrons primarily for aerobic respiration.",
           PUTATIVE),
        # [NiFe] Group 3: cofactor-coupled bidirectional
        _s("NiFe", "Group_3a", "F420-coupled",
           "Couples oxidation of H2 to reduction of F420 during methanogenesis. "
           "Physiologically reversible. [NiFeSe] variants.",
           CHARACTERIZED, _BIDIR, "f420_reducing"),
        _s("NiFe", "Group_3b", "NADP-coupled",
           "Couples oxidation of NADPH to evolution of H2. Physiologically reversible. "
           "May have sulfhydrogenase activity.",
           CHARACTERIZED, _BIDIR, "nadp_coupled"),
        _s("NiFe", "Group_3c", "Heterodisulfide reductase-linked",
           "Bifurcates electrons from H2 to heterodisulfide and Fdox in methanogens. [NiFeSe] variants.",
           CHARACTERIZED, _BIDIR, "heterodisulfide_reductase_linked", "bifurcating_hydrogenase"),
        _s("NiFe", "Group_3d", "NAD-coupled",
           "Interconverts electrons between H2 and NAD depending on cellular redox state.",
           CHARACTERIZED, _BIDIR, "nad_coupled"),
        # [NiFe] Group 4: respiratory H2-evolving
        _s("NiFe", "Group_4a", "Formate hydrogenlyase",
           "Couples formate oxidation to fermentative H2 evolution. May be H+-translocating.",
           CHARACTERIZED, _EVOLVING, "formate_coupled"),
        _s("NiFe", "Group_4b", "Formate-respiring",
           "Respires formate or carbon monoxide using H+ as electron acceptor. Na+-translocating via Mrp.",
           CHARACTERIZED, _EVOLVING, _ION, "formate_coupled"),
        _s("NiFe", "Group_4c", "Carbon monoxide-respiring",
           "Respires carbon monoxide using H+ as electron acceptor. H+-translocating.",
           CHARACTERIZED, _EVOLVING, _ION, "co_coupled"),
        _s("NiFe", "Group_4d", "Ferredoxin-coupled, Mrp-linked",
           "Couples Fdred oxidation to H+ reduction. Na+-translocating via Mrp complex.",
           CHARACTERIZED, _EVOLVING, _ION, _FD),
        _s("NiFe", "Group_4e", "Ferredoxin-coupled, Ech-type",
           "Couples Fdred oxidation to H+ reduction. Physiologically reversible via H+/Na+ translocation.",
           CHARACTERIZED, _EVOLVING, _ION, _FD, "ech_hydrogenase"),
        _s("NiFe", "Group_4f", "Formate-coupled (putative)",
           "Undetermined role. May couple formate oxidation to H2 evolution and H+ translocation.",
           PUTATIVE),
        _s("NiFe", "Group_4g", "Ferredoxin-coupled (putative)",
           "Undetermined role. May couple Fdred oxidation to proton reduction and H+/Na+ translocation.",
           PUTATIVE),
        _s("NiFe", "Group_4h", "Ferredoxin-coupled, Eha-type",
           "Couples Fdred oxidation to H+ reduction in anaplerotic processes. H+/Na+-translocating.",
           CHARACTERIZED, _EVOLVING, _ION, _FD),
        _s("NiFe", "Group_4i", "Ferredoxin-coupled, Ehb-type",
           "Couples Fdred oxidation to H+ reduction in anabolic processes. H+/Na+-translocating.",
           CHARACTERIZED, _EVOLVING, _ION, _FD),
        # [FeFe]: Group A subtypes need gene organization; sequence match stays structural
        _s("FeFe", "Group_A1", "Prototypical",
           "Couples ferredoxin oxidation to fermentative or photobiological H2 evolution.",
           CHARACTERIZED),
        _s("FeFe", "Group_A2", "Glutamate synthase-linked (putative)",
           "Undetermined role. May couple H2 oxidation to NAD reduction, generating reductant "
           "for glutamate synthase.",
           PUTATIVE),
        _s("FeFe", "Group_A3", "Bifurcating",
           "Reversibly bifurcates electrons from H2 to NAD and Fdox in anaerobic bacteria.",
           CHARACTERIZED),
        _s("FeFe", "Group_A4", "Formate dehydrogenase-linked",
           "Couples formate oxidation to H2 evolution. Some bifurcate electrons from H2 to "
           "ferredoxin and NADP.",
           CHARACTERIZED),
        _s("FeFe", "Group_B", "Colonic-type (putative)",
           "Undetermined role. May couple Fdred oxidation to fermentative H2 evolution.",
           PUTATIVE),
        _s("FeFe", "Group_C1", "Histidine kinase-linked (putative)",
           "Undetermined role. May sense H2 and regulate processes via histidine kinases.",
           PUTATIVE),
        _s("FeFe", "Group_C2", "Chemotactic (putative)",
           "Undetermined role. May sense H2 and regulate processes via methyl-accepting "
           "chemotaxis proteins.",
           PUTATIVE),
        _s("FeFe", "Group_C3", "Phosphatase-linked (putative)",
           "Undetermined role. May sense H2 and regulate processes via serine/threonine phosphatases.",
           PUTATIVE),
        # [Fe]
        _s("Fe", "Fe_only", "Methenyl-H4MPT dehydrogenase",
           "Reversibly couples H2 oxidation to 5,10-methenyltetrahydromethanopterin reduction.",
           CHARACTERIZED, "methanogen_hydrogenase"),
    ]
}

# Functional traits that sequence placement alone establishes only with
# genomic context (FeFe Group A subtypes). Recorded for documentation/tests.
CONTEXT_DEPENDENT_SUBGROUPS = frozenset(
    {("FeFe", "Group_A1"), ("FeFe", "Group_A2"), ("FeFe", "Group_A3"), ("FeFe", "Group_A4")}
)

# Terms this classifier emitted before the Table 1 audit; the refresh withdraws
# hydrogenase-derived copies of them.
RETIRED_TERMS = frozenset(
    {"nadp_reducing", "methyl_viologen_reducing", "monomeric_fefe", "sensory_hydrogenase"}
)


_REFERENCE_LABEL = re.compile(r"^\[(NiFe|FeFe)\]_(Group_\w+)$")


def parse_label(label: str) -> tuple[str | None, str | None]:
    """``[NiFe]_Group_4a`` -> ("NiFe", "Group_4a"); ``[Fe]`` -> ("Fe", "Fe_only")."""
    if label == "[Fe]":
        return "Fe", "Fe_only"
    match = _REFERENCE_LABEL.match(label)
    return (match.group(1), match.group(2)) if match else (None, None)


def group_of(hyd_type: str, subgroup: str) -> str:
    """HydDB group containing a subgroup: ``Group_4a`` -> ``Group_4``, ``Group_A3`` -> ``Group_A``."""
    if hyd_type == "Fe":
        return subgroup
    return subgroup[: len("Group_") + 1]


def lookup(hyd_type: str, subgroup: str) -> Subgroup | None:
    """Return the interpretation for an exact reference label, or None."""
    return SUBGROUPS.get((hyd_type, subgroup))
