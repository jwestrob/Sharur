# Embeddings and similarity search

Sharur stores one vector per protein and indexes them with FAISS for nearest-neighbor search. Similarity search only needs vectors, so any protein language model works.

## Choosing a model at ingest

Stage 06 of `sharur-ingest` embeds every protein with a Hugging Face protein encoder and mean-pools the residue vectors into one vector per protein. The default is ESM-2 8M, which is small enough for a laptop; choose another with `--embedding-model`:

```bash
sharur-ingest --input-dir genomes/ --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --embedding-model facebook/esm2_t33_650M_UR50D
```

The built-in stage loads models with `transformers.AutoModel` and `AutoTokenizer`, which covers the ESM-2 family and other encoders published in that form. Proteins longer than the model's context window are truncated; the embedding manifest counts them.

## Bringing your own embeddings

For models with their own inference code (ESM-C, ProtT5, structure-aware models, or anything else), compute the vectors elsewhere and write one HDF5 file:

| Dataset or attribute | Contents |
|---|---|
| `protein_ids` | 1-D UTF-8 strings matching `proteins.protein_id` in the database |
| `embeddings` | 2-D float array, one row per protein, any dimension |
| `model_name` (attribute, optional) | Recorded in the index manifest |

Place it at `data/my_dataset/embeddings/protein_embeddings.h5` and build the index:

```bash
sharur build-vector-index --embeddings data/my_dataset/embeddings/protein_embeddings.h5
```

`find_similar` and the other similarity operators then use your vectors. Queries must come from the same model as the index; the index records its dimension and rejects mismatched queries.

## Synteny

ELSA builds conserved gene blocks from the same per-protein embeddings and writes them to a run-scoped `synteny.duckdb` beside the dataset. See [Tools](../tools_reference.md).
