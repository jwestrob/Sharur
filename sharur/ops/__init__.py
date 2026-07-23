"""Sharur Ops — agent coordination layer.

Use OpsStore for local, single-process SQLite access:
    from sharur.ops.store import OpsStore
    ops = OpsStore("sharur_ops.db", agent_id="my_agent")

Use SharurOps for distributed access through one database-owning server:
    from sharur.ops.client import SharurOps
    ops = SharurOps("http://localhost:8811", agent_id="my_agent")

Both expose the same core record/task API. HTTP mode adds credential-derived
roles, replayable events, bounded connection pooling, and administrative
endpoints.
"""
