# Sharur Quick Reference Card

## Critical Patterns

### 1. Direct Store Queries

```python
# Use DuckDB store directly for protein queries
store = session.db
rows = store.execute(
    "SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence "
    "FROM proteins p "
    "JOIN annotations a ON p.protein_id = a.protein_id "
    "WHERE a.source = 'pfam' AND a.accession = ? "
    "AND EXISTS (SELECT 1 FROM contigs c WHERE c.contig_id = p.contig_id AND c.taxonomy LIKE ?)",
    ['PF00142', '%Archaea%'],
)
```

### 2. Strand-Aware Coordinates

```python
def get_upstream(protein, n_bp):
    if protein.strand == "+":
        return (protein.start - n_bp, protein.start)
    else:  # "-" strand
        return (protein.end, protein.end + n_bp)  # HIGHER coordinates!
```

### 3. Circular Contig Handling

```python
def normalize(coord, length):
    return coord % length

def window_wraps(start, end, length):
    """True if window crosses origin on circular contig"""
    return start > end

# Query: WHERE (start <= pos AND pos < end) OR 
#              (wraps AND (pos >= start OR pos < end))
```

### 4. Focus Stack Resolution

```python
RESOLUTION_MAP = {
    "it": lambda stack: next((e for e in reversed(stack) if e.type == "protein"), None),
    "them": lambda stack: next((e for e in reversed(stack) if e.type == "set"), None),
    "those": lambda stack: next((e for e in reversed(stack) if e.type == "set"), None),
    "that locus": lambda stack: next((e for e in reversed(stack) if e.type == "locus"), None),
}
```

### 5. Tool Result Format

```python
@dataclass
class ToolResult:
    success: bool        # Did execution succeed?
    data: Any            # The actual results
    summary: str         # "Found 47 proteins with PF00142"
    count: int           # len(data) if applicable
    truncated: bool      # True if limit was hit
    warnings: list[str]  # Reliability concerns
    suggestions: list[str]  # Follow-up ideas
```

### 6. Session State Updates

Every tool execution should:
```python
def execute(self, params, session):
    results = self._do_search(params)
    
    # Update focus stack
    if results:
        session.push_focus("protein_set", results.to_id_list(), f"{len(results)} proteins")
    
    # Log provenance
    session.log_query(
        query=params.original_query,
        tool_calls=[{"tool": self.name, "params": params.dict()}],
        results_summary=results.summary,
        duration_ms=elapsed
    )
    
    return ToolResult(success=True, data=results, ...)
```

## Key Tables

### Loci Table (CRISPR Arrays, BGCs, etc.)
```sql
-- CRISPR arrays are detected by MinCED (stage 05c) and stored in loci table
SELECT * FROM loci WHERE locus_type = 'crispr_array';

-- BGCs from GECCO (stage 05a)
SELECT * FROM loci WHERE locus_type = 'bgc';

-- Schema:
--   locus_id     VARCHAR PRIMARY KEY
--   locus_type   VARCHAR  -- 'crispr_array', 'bgc', 'prophage'
--   contig_id    VARCHAR REFERENCES contigs
--   start        INTEGER
--   end_coord    INTEGER
--   confidence   FLOAT
--   metadata     JSON  -- Contains repeat sequences, etc.
```

**Note:** `visualize_neighborhood()` automatically overlays CRISPR arrays when present.

---

## SQL Patterns

### Spatial Window Query
```sql
SELECT p2.* FROM proteins p1
JOIN proteins p2 ON p1.contig_id = p2.contig_id
JOIN contigs c ON p1.contig_id = c.contig_id
WHERE p1.protein_id = :anchor
AND (
    -- Linear contig
    (NOT c.is_circular AND p2.start >= p1.start - :window AND p2.end <= p1.end + :window)
    OR
    -- Circular contig (wrap-around case)
    (c.is_circular AND (
        (p2.start >= :wrapped_start OR p2.end <= :wrapped_end)
    ))
)
ORDER BY p2.start;
```

### Domain Search with Taxonomy
```sql
SELECT p.* FROM proteins p
JOIN annotations a ON p.protein_id = a.protein_id
JOIN bins b ON p.bin_id = b.bin_id
WHERE a.source = 'pfam' AND a.accession = :domain
AND b.taxonomy LIKE :tax_pattern || '%';
```

### Operon Expansion
```sql
WITH RECURSIVE operon AS (
    -- Seed: the starting gene
    SELECT protein_id, contig_id, start, end, strand, 0 as depth
    FROM proteins WHERE protein_id = :seed
    
    UNION ALL
    
    -- Expand: adjacent genes on same strand within gap threshold
    SELECT p.protein_id, p.contig_id, p.start, p.end, p.strand, o.depth + 1
    FROM proteins p
    JOIN operon o ON p.contig_id = o.contig_id AND p.strand = o.strand
    WHERE (p.start - o.end BETWEEN 0 AND :gap_threshold 
           OR o.start - p.end BETWEEN 0 AND :gap_threshold)
    AND p.protein_id != o.protein_id
    AND o.depth < 50  -- Safety limit
)
SELECT DISTINCT * FROM operon;
```

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| "No results" for obvious query | Taxonomy mismatch | Check GTDB format vs query |
| Wrong neighbors | Forgot strand | Use strand-aware upstream/downstream |
| Infinite loop | Circular contig | Normalize coordinates |
| Empty focus | Didn't push results | Always push after tool returns |
## File Dependencies

```
types.py          <- Everything imports from here
    ↓
storage/*.py      <- Needs types
    ↓
operators/*.py    <- Needs storage + types
    ↓
predicates/*.py   <- Needs storage + types
    ↓
session.py        <- Needs types + operators
```

## Testing Checklist

- [ ] Circular contig: window crosses origin
- [ ] Minus strand: upstream is higher coords
- [ ] Empty results: graceful handling
- [ ] Focus resolution: "it" -> most recent protein
- [ ] Working sets: CRUD operations
- [ ] Long query: timeout handling
- [ ] SQL injection: parameterized queries
- [ ] Missing vector store: graceful degradation
