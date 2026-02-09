"""
Predicate mappings from annotation sources.

Maps PFAM, KEGG, CAZy, and VOGdb annotations to predicates.
"""

from sharur.predicates.mappings.pfam_map import PFAM_TO_PREDICATES, PFAM_PATTERNS
from sharur.predicates.mappings.kegg_map import KEGG_TO_PREDICATES, EC_TO_PREDICATES
from sharur.predicates.mappings.cazy_map import CAZY_TO_PREDICATES, CAZY_FAMILY_PATTERNS
from sharur.predicates.mappings.vog_map import (
    VOG_CATEGORY_TO_PREDICATES,
    VOG_DESCRIPTION_PATTERNS,
    VOG_DIRECT_MAPPINGS,
)

__all__ = [
    "PFAM_TO_PREDICATES",
    "PFAM_PATTERNS",
    "KEGG_TO_PREDICATES",
    "EC_TO_PREDICATES",
    "CAZY_TO_PREDICATES",
    "CAZY_FAMILY_PATTERNS",
    "VOG_CATEGORY_TO_PREDICATES",
    "VOG_DESCRIPTION_PATTERNS",
    "VOG_DIRECT_MAPPINGS",
]
