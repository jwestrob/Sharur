"""HydDB-based hydrogenase subgroup assignment, interpretation, and refresh."""

from sharur.hydrogenase.classifier import (
    CLASSIFIER_VERSION,
    Classification,
    HydrogenaseSearchError,
    classify,
    classify_database,
    write_classifications,
)
from sharur.hydrogenase.refresh import RefreshValidationError, refresh_hydrogenases
from sharur.hydrogenase.subgroups import SUBGROUPS, Subgroup, lookup


__all__ = [
    "CLASSIFIER_VERSION",
    "SUBGROUPS",
    "Classification",
    "HydrogenaseSearchError",
    "RefreshValidationError",
    "Subgroup",
    "classify",
    "classify_database",
    "lookup",
    "refresh_hydrogenases",
    "write_classifications",
]
