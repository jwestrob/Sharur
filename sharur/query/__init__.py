"""Read-only, bounded analytical data plane for coordinated Sharur agents."""

from sharur.query.client import SharurQuery
from sharur.query.staging import StagedDatabase, stage_database


__all__ = ["SharurQuery", "StagedDatabase", "stage_database"]
