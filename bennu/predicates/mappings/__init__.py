"""
Predicate mappings from annotation sources.

Maps PFAM, KEGG, CAZy, and VOGdb annotations to predicates.
"""

from bennu.predicates.mappings.pfam_map import PFAM_TO_PREDICATES, PFAM_PATTERNS
from bennu.predicates.mappings.kegg_map import KEGG_TO_PREDICATES, EC_TO_PREDICATES
from bennu.predicates.mappings.cazy_map import CAZY_TO_PREDICATES, CAZY_FAMILY_PATTERNS
from bennu.predicates.mappings.vog_map import (
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
