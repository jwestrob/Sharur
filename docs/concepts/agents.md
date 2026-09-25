# Working with agents

Sharur's interface is a command line and a Python API whose outputs are bounded, typed and free of raw sequences. Any agent that can run shell commands or Python can use it.

| Agent | Setup |
|---|---|
| **Claude Code** | Open the repository. Claude Code reads `CLAUDE.md` and the skills in `.claude/skills/`, invoked as `/survey`, `/explore`, `/characterize`, `/hydrogenase` and so on. |
| **Codex** | Open the repository. Codex reads `AGENTS.md`, which is the same file as `CLAUDE.md`. |
| **Other agents, including open-weight models** | Give the agent shell or Python access and point it at `CLAUDE.md` for the rules and at `sharur describe`, `card` and `why` for bounded summaries. |

## The rules agents follow

`CLAUDE.md` holds the rules every agent works under. The main ones:

- **Observed versus named.** Report the domains a protein carries as observations; name a system, family or pathway only when a curated caller made that call.
- **Every number has a query.** Each count or identifier in a finding carries the query that reproduces it.
- **No sequences in model-visible text.** Agents refer to proteins by identifier; sequence work happens in tools.
- **Incomplete genomes.** "Not detected in this MAG" instead of "the genome lacks".

## Skills

| Skill | Use |
|---|---|
| `survey` | Systematic survey of a dataset |
| `explore` | Open-ended discovery around loci of interest |
| `characterize` | One protein or locus in depth, including structure |
| `metabolism`, `pathway` | Pathway reconstruction and module completeness |
| `defense`, `prophage` | Defense systems and viral elements |
| `hydrogenase` | NiFe and FeFe hydrogenase validation |
| `compare` | Cross-genome comparison |
| `synteny` | ELSA conserved gene blocks |
| `foldseek` | Structure prediction and structural search |
| `literature` | Literature lookups for ambiguous functions |
| `reviewer_2` | Adversarial check of claims |
| `coordinator`, `atlas`, `brainstorm`, `query`, `visualize` | Orchestration, genome-by-genome reading, synthesis, quick queries, figures |

## Campaigns

For reading every genome in a large dataset, Sharur coordinates many agents through three services: `sharur-ops` (tasks and findings), `sharur-query` (one shared read-only database) and `sharur-review` (layered review of candidates). `sharur-worker` runs headless Claude Code (`claude -p`) or Codex (`codex exec`) sessions against those tasks. See [Query service](../query_service.md), [Agent coordination](../agent_ops_spec.md) and [Review workflow](../review_workflow.md).
