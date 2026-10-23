"""DefenseFinder per-hit selection with MacSyFinder/DefenseFinder default thresholds.

DefenseFinder calls systems by running MacSyFinder over the defense-finder-models
package. MacSyFinder decides which HMM domain hits enter system calling in two steps,
and this module reproduces both on a permissive hits table (one row per domain, as
written by Aksha/Astra ``search``):

1. Search threshold. Profiles that carry a GA line are searched with ``hmmsearch
   --cut_ga``: the full-sequence score must reach GA1 and the domain score GA2.
   Profiles without GA are searched with ``-E e_value_search`` (0.1).
2. Hit selection. Every reported domain is then kept only when its independent
   E-value is at most ``i_evalue_sel`` (1e-3) and its profile coverage,
   ``(hmm_to - hmm_from + 1) / LENG``, is at least ``coverage_profile``.

Defaults, read from the installed packages (2026-10-23), never from memory:

* macsylib 1.0.4 (the library under MacSyFinder 2.1.4), ``macsylib/config.py``
  ``MacsyDefaults``: ``cut_ga=True``, ``e_value_search=0.1``, ``i_evalue_sel=0.001``,
  ``coverage_profile=0.5``.
* ``macsylib/profile.py`` ``Profile.execute``: ``--cut_ga`` when the profile has a GA
  line, otherwise ``-E e_value_search``.
* ``macsylib/report.py`` ``HMMReport._parse_hmm_body``: keeps a domain only when
  ``i_evalue <= i_evalue_sel`` and ``profile coverage >= coverage_profile``. The GA
  search threshold and this selection apply together.
* mdmparis-defense-finder 2.0.1, ``defense_finder_cli/main.py`` ``run``: ``--coverage``
  defaults to 0.4 and reaches MacSyFinder as ``--coverage-profile 0.4``;
  ``defense_finder/__init__.py`` leaves ``cut_ga`` at MacSyFinder's default.

The default coverage here is therefore DefenseFinder's 0.4; MacSyFinder's own
default (0.5) is kept as a named constant for runs that want it.

Profile lengths and GA values come from the installed model set's HMM files. Hits
name their profile by the HMM ``NAME`` field, which differs from the file stem for
about a quarter of defense-finder-models profiles, so profiles are indexed by both.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd


if TYPE_CHECKING:
    import duckdb


logger = logging.getLogger(__name__)

SOURCE = "defensefinder"
MODEL_FAMILY = "defense-finder-models"

# MacSyFinder defaults (macsylib 1.0.4 config.py MacsyDefaults).
MACSYFINDER_CUT_GA = True
MACSYFINDER_E_VALUE_SEARCH = 0.1
MACSYFINDER_I_EVALUE_SEL = 1e-3
MACSYFINDER_COVERAGE_PROFILE = 0.5
# DefenseFinder 2.0.1 passes --coverage-profile 0.4 (defense_finder_cli/main.py).
DEFENSEFINDER_COVERAGE_PROFILE = 0.4

DEFAULT_I_EVALUE_SEL = MACSYFINDER_I_EVALUE_SEL
DEFAULT_COVERAGE_PROFILE = DEFENSEFINDER_COVERAGE_PROFILE

SELECTION_RULE = (
    "MacSyFinder/DefenseFinder: GA1/GA2 bit thresholds where the profile has GA "
    "(full-sequence E <= e_value_search otherwise), then i-Evalue <= i_evalue_sel "
    "and profile coverage >= coverage_profile"
)

# Sub-families DefenseFinder 2.0.1 runs by default (defense_finder/__init__.py run()):
# DefenseFinder and RM. Anti-defense (ADF) models and their profiles run only with -a.
DEFAULT_SUBFAMILIES: tuple[str, ...] = ("DefenseFinder", "RM")

#: Environment variable naming a defense-finder-models directory (or its profiles/).
MODELS_ENV = "SHARUR_DEFENSEFINDER_MODELS"

_REQUIRED_COLUMNS = (
    "sequence_id", "hmm_name", "bitscore", "evalue", "i_evalue", "dom_bitscore",
    "hmm_from", "hmm_to",
)

# DefenseFinder accessions with numbered FAM variants: trailing _FAM_N or _FAMN.
#   RM_Type_II__Type_II_MTases_FAM_1  -> group RM_Type_II__Type_II_MTases
#   RM_Type_III__Type_III_MTases_FAM10 -> group RM_Type_III__Type_III_MTases
# Component suffixes such as MkoA_B carry no leading _FAM and stay distinct.
FAM_VARIANT_RE = re.compile(r"_FAM_?\d+$")


@dataclass(frozen=True)
class ProfileInfo:
    """Header fields of one HMM profile."""

    name: str
    stem: str
    length: int
    ga_seq: float | None = None
    ga_dom: float | None = None


@dataclass
class ModelSet:
    """An installed defense-finder-models directory (or a flat profile directory)."""

    profiles_dir: Path
    root: Path | None
    version: str | None
    profiles: list[ProfileInfo]

    @property
    def models_dir(self) -> Path | None:
        """MacSyFinder models root (parent of ``defense-finder-models``) for co-location."""
        if self.root is not None and (self.root / "definitions").is_dir():
            return self.root.parent
        return None

    @property
    def family(self) -> str:
        return self.root.name if self.root is not None else MODEL_FAMILY

    def frame(self) -> pd.DataFrame:
        """Profile lookup keyed by HMM NAME, with file stems as a fallback key."""
        by_key: dict[str, ProfileInfo] = {}
        for info in self.profiles:
            by_key.setdefault(info.name, info)
        for info in self.profiles:
            by_key.setdefault(info.stem, info)
        return pd.DataFrame(
            {
                "hmm_name": list(by_key),
                "_profile_length": [p.length for p in by_key.values()],
                "_ga_seq": [np.nan if p.ga_seq is None else p.ga_seq for p in by_key.values()],
                "_ga_dom": [np.nan if p.ga_dom is None else p.ga_dom for p in by_key.values()],
            }
        )


@dataclass
class SelectionParams:
    use_ga: bool = MACSYFINDER_CUT_GA
    e_value_search: float = MACSYFINDER_E_VALUE_SEARCH
    i_evalue_sel: float = DEFAULT_I_EVALUE_SEL
    coverage_profile: float = DEFAULT_COVERAGE_PROFILE
    # Load hits only for profiles that these definitions/ sub-families reference: the
    # models DefenseFinder runs by default, which the co-location caller also runs.
    subfamilies: tuple[str, ...] = DEFAULT_SUBFAMILIES


@dataclass
class SelectionStats:
    raw_rows: int = 0
    unknown_profile: int = 0
    out_of_scope_profile: int = 0
    missing_alignment: int = 0
    failed_threshold: int = 0
    failed_i_evalue: int = 0
    failed_coverage: int = 0
    selected_rows: int = 0
    fam_dedup_removed: int = 0
    loaded_rows: int = 0
    unknown_profile_names: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def parse_hmm_header(path: Path) -> ProfileInfo | None:
    """Read NAME, LENG and GA from one HMM file header (first model only)."""
    name = length = None
    ga_seq = ga_dom = None
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("NAME"):
                parts = line.split()
                name = parts[1] if len(parts) > 1 else None
            elif line.startswith("LENG"):
                length = int(line.split()[1])
            elif line.startswith("GA"):
                parts = line.replace(";", " ").split()
                try:
                    ga_seq, ga_dom = float(parts[1]), float(parts[2])
                except (IndexError, ValueError):
                    # MacSyFinder ignores a malformed GA line (profile.py _profile_features).
                    ga_seq = ga_dom = None
            elif line.startswith(("STATS", "HMM ")):
                break
    if name is None or length is None:
        return None
    return ProfileInfo(name=name, stem=path.stem, length=length, ga_seq=ga_seq, ga_dom=ga_dom)


def read_model_version(root: Path | None) -> str | None:
    """``vers`` from a defense-finder-models ``metadata.yml``; None when absent."""
    if root is None:
        return None
    metadata = root / "metadata.yml"
    if not metadata.is_file():
        return None
    for line in metadata.read_text(errors="replace").splitlines():
        match = re.match(r"^vers\s*:\s*['\"]?([^'\"\s]+)", line)
        if match:
            return match.group(1)
    return None


def load_model_set(path: str | Path) -> ModelSet:
    """Load a defense-finder-models directory, its ``profiles/``, or a flat profile directory."""
    path = Path(path).expanduser().resolve()
    if (path / "profiles").is_dir():
        root, profiles_dir = path, path / "profiles"
    elif path.is_dir() and path.name == "profiles" and (path.parent / "metadata.yml").is_file():
        root, profiles_dir = path.parent, path
    elif path.is_dir():
        root, profiles_dir = None, path
    else:
        raise FileNotFoundError(f"DefenseFinder model directory not found: {path}")
    profiles = [
        info for hmm in sorted(profiles_dir.glob("*.hmm")) if (info := parse_hmm_header(hmm))
    ]
    if not profiles:
        raise FileNotFoundError(f"No readable .hmm profiles in {profiles_dir}")
    return ModelSet(profiles_dir=profiles_dir, root=root, version=read_model_version(root),
                    profiles=profiles)


def find_model_set(explicit: str | Path | None = None) -> ModelSet:
    """Resolve the installed model set: explicit path, $SHARUR_DEFENSEFINDER_MODELS, MacSyFinder home."""
    candidates: list[Path] = []
    if explicit is not None:
        return load_model_set(explicit)
    if os.environ.get(MODELS_ENV):
        return load_model_set(os.environ[MODELS_ENV])
    candidates = [
        Path.home() / ".macsyfinder" / "models" / MODEL_FAMILY,
        Path.home() / ".mdmlab" / "macsyfinder" / "models" / MODEL_FAMILY,
    ]
    for candidate in candidates:
        if (candidate / "profiles").is_dir():
            return load_model_set(candidate)
    raise FileNotFoundError(
        "DefenseFinder models not found; pass the defense-finder-models directory explicitly "
        f"or set {MODELS_ENV} (searched {', '.join(str(c) for c in candidates)})"
    )


_GENE_NAME_RE = re.compile(r'<gene\s+name="([^"]+)"')


def referenced_profiles(model_set: ModelSet, subfamilies: tuple[str, ...]) -> set[str] | None:
    """Profile NAMEs and stems the given sub-families' models reference (exchangeables too).

    None when the model set carries no definitions (a flat profile directory) or no
    sub-families are named: every profile is then in scope.
    """
    if not subfamilies or model_set.root is None:
        return None
    definitions = model_set.root / "definitions"
    if not definitions.is_dir():
        return None
    from sharur.colocation import _definition_roots

    stems: set[str] = set()
    for root in _definition_roots(definitions, tuple(subfamilies)):
        for xml in root.rglob("*.xml"):
            stems.update(_GENE_NAME_RE.findall(xml.read_text(errors="replace")))
    names = {p.name for p in model_set.profiles if p.stem in stems}
    return stems | names


def read_hits(path: str | Path) -> pd.DataFrame:
    """Read an Aksha/Astra per-domain hits table."""
    return pd.read_csv(
        path, sep="\t",
        dtype={"sequence_id": str, "hmm_name": str},
        na_values=[""], keep_default_na=False,
    )


def select_hits(
    hits: pd.DataFrame,
    model_set: ModelSet,
    params: SelectionParams | None = None,
) -> tuple[pd.DataFrame, SelectionStats, pd.DataFrame]:
    """Apply MacSyFinder/DefenseFinder hit selection to a per-domain hits table.

    Returns (selected rows with the input columns, stats, per-row gate table). The gate
    table carries ``profile_coverage`` and boolean ``pass_*`` columns for every input
    row, for sensitivity summaries.
    """
    params = params or SelectionParams()
    missing = [c for c in _REQUIRED_COLUMNS if c not in hits.columns]
    if missing:
        raise ValueError(f"DefenseFinder hits table lacks required columns: {missing}")
    stats = SelectionStats(raw_rows=len(hits))
    if hits.empty:
        return hits.copy(), stats, pd.DataFrame(index=hits.index)

    joined = hits.merge(model_set.frame(), on="hmm_name", how="left", validate="many_to_one")
    joined.index = hits.index
    known = joined["_profile_length"].notna()
    scope = referenced_profiles(model_set, params.subfamilies)
    in_scope = joined["hmm_name"].isin(scope) if scope is not None else known | True
    considered = known & in_scope
    hmm_from = pd.to_numeric(joined["hmm_from"], errors="coerce")
    hmm_to = pd.to_numeric(joined["hmm_to"], errors="coerce")
    coverage = (hmm_to - hmm_from + 1) / joined["_profile_length"]
    has_alignment = coverage.notna()
    bitscore = pd.to_numeric(joined["bitscore"], errors="coerce")
    dom_bitscore = pd.to_numeric(joined["dom_bitscore"], errors="coerce")
    full_evalue = pd.to_numeric(joined["evalue"], errors="coerce")
    i_evalue = pd.to_numeric(joined["i_evalue"], errors="coerce")

    has_ga = joined["_ga_seq"].notna() & joined["_ga_dom"].notna()
    if params.use_ga:
        ga_pass = (bitscore >= joined["_ga_seq"]) & (dom_bitscore >= joined["_ga_dom"])
        threshold_pass = (has_ga & ga_pass) | (~has_ga & (full_evalue <= params.e_value_search))
    else:
        threshold_pass = full_evalue <= params.e_value_search
    ie_pass = i_evalue <= params.i_evalue_sel
    cov_pass = coverage >= params.coverage_profile
    keep = considered & has_alignment & threshold_pass & ie_pass & cov_pass

    stats.unknown_profile = int((~known).sum())
    stats.out_of_scope_profile = int((known & ~in_scope).sum())
    stats.missing_alignment = int((considered & ~has_alignment).sum())
    stats.failed_threshold = int((considered & has_alignment & ~threshold_pass).sum())
    stats.failed_i_evalue = int((considered & has_alignment & threshold_pass & ~ie_pass).sum())
    stats.failed_coverage = int(
        (considered & has_alignment & threshold_pass & ie_pass & ~cov_pass).sum()
    )
    stats.selected_rows = int(keep.sum())
    stats.unknown_profile_names = sorted(joined.loc[~known, "hmm_name"].dropna().unique())[:20]

    gates = pd.DataFrame(
        {
            "hmm_name": joined["hmm_name"],
            # known profile referenced by the called sub-families
            "known_profile": considered,
            "in_scope": in_scope,
            "has_ga": has_ga,
            "profile_coverage": coverage,
            "pass_threshold": threshold_pass,
            "pass_ga_only": considered & has_ga & ga_pass if params.use_ga else considered & False,
            "pass_i_evalue": ie_pass,
            "pass_coverage": cov_pass,
            "selected": keep,
        },
        index=hits.index,
    )
    return hits.loc[keep].copy(), stats, gates


def dedup_fam_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Keep the best hit per (protein, FAM-variant group) for numbered FAM profiles.

    DefenseFinder has ~40 ``Type_II_MTases_FAM_*`` profiles that all detect variants of
    one Rossmann methyltransferase fold, so one protein can match many of them. Rows are
    grouped by ``protein_id`` and the accession with a trailing ``_FAM_N`` removed; the
    highest-score row (then lowest e-value) is kept. Accessions without the suffix pass
    through unchanged.

    Args:
        df: annotation rows with at least ``protein_id``, ``accession``, ``score``, ``evalue``.
    """
    if df.empty:
        return df
    df = df.copy()
    df["_dedup_group"] = df["accession"].apply(
        lambda x: FAM_VARIANT_RE.sub("", x) if isinstance(x, str) else x
    )
    has_fam = df["_dedup_group"] != df["accession"]
    if not has_fam.any():
        return df.drop(columns=["_dedup_group"])

    no_fam_df = df[~has_fam].drop(columns=["_dedup_group"])
    fam_df = df[has_fam].sort_values(["score", "evalue"], ascending=[False, True],
                                     na_position="last")
    fam_deduped = fam_df.drop_duplicates(subset=["protein_id", "_dedup_group"],
                                         keep="first").drop(columns=["_dedup_group"])
    n_removed = len(fam_df) - len(fam_deduped)
    if n_removed:
        logger.info(
            "  defensefinder: deduplicated %s redundant FAM variant hits (%s -> %s FAM rows)",
            f"{n_removed:,}", f"{len(fam_df):,}", f"{len(fam_deduped):,}",
        )
    return pd.concat([no_fam_df, fam_deduped], ignore_index=True)


ANNOTATION_COLUMNS = [
    "protein_id", "source", "accession", "name", "description", "evalue", "score",
    "start_aa", "end_aa",
]


def annotation_rows(selected: pd.DataFrame) -> pd.DataFrame:
    """Selected hits as ``annotations`` rows, matching Stage 07's column conventions.

    ``evalue`` is the full-sequence E-value (filled from i/c-Evalue when absent), ``score``
    the domain bit score (the score MacSyFinder uses), coordinates the alignment span
    (filled from the envelope). ``name`` repeats the profile NAME; ``description`` is empty.
    """
    def column(name: str) -> pd.Series:
        if name in selected:
            return pd.to_numeric(selected[name], errors="coerce")
        return pd.Series(np.nan, index=selected.index)

    rows = pd.DataFrame(
        {
            "protein_id": selected["sequence_id"].astype(str),
            "source": SOURCE,
            "accession": selected["hmm_name"].astype(str),
            "name": selected["hmm_name"].astype(str),
            "description": "",
            "evalue": column("evalue").fillna(column("i_evalue")).fillna(column("c_evalue")),
            "score": column("dom_bitscore").fillna(column("bitscore")),
            "start_aa": column("ali_from").fillna(column("env_from")),
            "end_aa": column("ali_to").fillna(column("env_to")),
        },
        index=selected.index,
    )
    return rows[ANNOTATION_COLUMNS]


def build_annotations(
    hits: pd.DataFrame,
    model_set: ModelSet,
    params: SelectionParams | None = None,
) -> tuple[pd.DataFrame, SelectionStats]:
    """Selection, row conversion and FAM-variant dedup: the rows Stage 07 loads."""
    selected, stats, _ = select_hits(hits, model_set, params)
    rows = annotation_rows(selected)
    deduped = dedup_fam_variants(rows)
    stats.fam_dedup_removed = len(rows) - len(deduped)
    stats.loaded_rows = len(deduped)
    if stats.unknown_profile:
        logger.warning(
            "defensefinder: %s hit rows name %s profile(s) absent from %s (e.g. %s); "
            "they cannot be thresholded and were not loaded",
            f"{stats.unknown_profile:,}", len(stats.unknown_profile_names) or "some",
            model_set.profiles_dir, ", ".join(stats.unknown_profile_names[:5]),
        )
    return deduped.reset_index(drop=True), stats


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


SOURCE_VERSIONS_TABLE = "annotation_source_versions"
_SOURCE_VERSIONS_SQL = f"""
CREATE TABLE IF NOT EXISTS {SOURCE_VERSIONS_TABLE} (
    source VARCHAR PRIMARY KEY,
    model_set VARCHAR,
    model_version VARCHAR,
    profiles_path VARCHAR,
    profile_count INTEGER,
    selection_rule VARCHAR,
    parameters VARCHAR,
    hits_path VARCHAR,
    hits_sha256 VARCHAR,
    raw_rows BIGINT,
    selected_rows BIGINT,
    loaded_rows BIGINT,
    loaded_at TIMESTAMP
)
"""


def record_source_version(
    conn: duckdb.DuckDBPyConnection,
    model_set: ModelSet,
    params: SelectionParams,
    stats: SelectionStats,
    hits_path: str | Path | None,
    hits_sha256: str | None = None,
) -> None:
    """Record the model set, selection parameters and hit counts behind the loaded rows."""
    conn.execute(_SOURCE_VERSIONS_SQL)
    conn.execute(
        f"INSERT OR REPLACE INTO {SOURCE_VERSIONS_TABLE} VALUES "
        "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            SOURCE, model_set.family, model_set.version, str(model_set.profiles_dir),
            len(model_set.profiles), SELECTION_RULE, json.dumps(asdict(params), sort_keys=True),
            str(hits_path) if hits_path is not None else None, hits_sha256,
            stats.raw_rows, stats.selected_rows, stats.loaded_rows,
            datetime.now(timezone.utc).replace(tzinfo=None),
        ],
    )
