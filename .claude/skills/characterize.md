# Characterize Skill (Deprecated)

**This skill has been folded into the explore and foldseek workflows.**

- For protein characterization: use `/explore` with targeted hypothesis testing
- For structure prediction and Foldseek search: use `/foldseek`
- For literature research on hits: use `/literature`
- For genomic context: the explore agent handles this natively via `get_neighborhood()` and co-annotation queries

The original characterize pipeline (context -> structure -> foldseek -> literature -> hypothesis) is now the natural workflow of an explore agent dispatching specialist sub-agents as needed.
