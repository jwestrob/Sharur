"""Sharur Ops — agent coordination layer.

Use OpsStore for direct SQLite access (no server needed):
    from sharur.ops.store import OpsStore
    ops = OpsStore("sharur_ops.db", agent_id="my_agent")

Use SharurOps for HTTP access (requires server running):
    from sharur.ops.client import SharurOps
    ops = SharurOps("http://localhost:8811", agent_id="my_agent")

Both expose the same API: finding(), recent_findings(), hypothesis(), etc.
"""
