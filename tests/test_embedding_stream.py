"""Memory-bounded Stage 06 FASTA and HDF5 stream tests."""

import h5py
import numpy as np
import pytest

from sharur.ingest.embedding_stream import (
    iter_stage03_proteins,
    stream_embeddings_to_h5,
)


def test_stage03_proteins_stream_in_deterministic_genome_order(tmp_path):
    genomes = tmp_path / "stage03" / "genomes"
    (genomes / "genome_b").mkdir(parents=True)
    (genomes / "genome_a").mkdir()
    (genomes / "all_protein_symlinks").mkdir()
    (genomes / "genome_b" / "b.faa").write_text(">b1 description\nMKT*\n")
    (genomes / "genome_a" / "a.faa").write_text(">a1\nMA\n>a2\nVVV*\n")
    (genomes / "all_protein_symlinks" / "ignored.faa").write_text(">ignored\nM\n")

    records = list(iter_stage03_proteins(tmp_path / "stage03"))

    assert records == [
        ("a1", "MA", "genome_a"),
        ("a2", "VVV", "genome_a"),
        ("b1", "MKT", "genome_b"),
    ]


def test_embeddings_stream_extends_h5_and_accumulates_statistics(tmp_path):
    records = (
        record
        for record in [
            ("p1", "MA", "g1"),
            ("p2", "VVV", "g1"),
            ("p3", "LLLL", "g2"),
        ]
    )
    batch_sizes = []

    def embed(sequences):
        batch_sizes.append(len(sequences))
        return np.asarray([[len(sequence), 1.0] for sequence in sequences], dtype=np.float32)

    h5_path = tmp_path / "embeddings.h5"
    stats = stream_embeddings_to_h5(
        records,
        h5_path,
        embedding_dim=2,
        model_name="test-model",
        batch_size=2,
        embed_batch=embed,
    )

    assert batch_sizes == [2, 1]
    assert stats.to_dict() == {
        "total_proteins": 3,
        "genomes_processed": 2,
        "min_protein_length": 2,
        "max_protein_length": 4,
        "mean_protein_length": 3.0,
        "truncated_sequences": 0,
        "genome_protein_counts": {"g1": 2, "g2": 1},
    }
    with h5py.File(h5_path, "r") as handle:
        assert [value.decode() for value in handle["protein_ids"][:]] == ["p1", "p2", "p3"]
        np.testing.assert_allclose(
            handle["embeddings"][:],
            np.asarray([[2, 1], [3, 1], [4, 1]], dtype=np.float32),
        )
        assert handle.attrs["num_proteins"] == 3
        assert handle.attrs["model_name"] == "test-model"


def test_embeddings_stream_counts_sequences_subject_to_truncation(tmp_path):
    stats = stream_embeddings_to_h5(
        [("short", "MA", "g1"), ("long", "VVVV", "g1")],
        tmp_path / "embeddings.h5",
        embedding_dim=2,
        model_name="test-model",
        batch_size=2,
        truncation_residue_limit=3,
        embed_batch=lambda sequences: np.ones((len(sequences), 2), dtype=np.float32),
    )

    assert stats.truncated_sequences == 1


def test_embeddings_stream_rejects_empty_sequences(tmp_path):
    with pytest.raises(ValueError, match="empty sequence"):
        stream_embeddings_to_h5(
            [("empty", "", "g1")],
            tmp_path / "empty.h5",
            embedding_dim=2,
            model_name="test-model",
            batch_size=1,
            embed_batch=lambda sequences: np.ones((len(sequences), 2), dtype=np.float32),
        )


@pytest.mark.parametrize(
    ("bad_batch", "message"),
    [
        (np.zeros((1, 2), dtype=np.float32), "zero vector"),
        (np.asarray([[np.nan, 1.0]], dtype=np.float32), "NaN or infinite"),
        (np.ones((1, 3), dtype=np.float32), "shape"),
    ],
)
def test_embeddings_stream_rejects_invalid_model_output(tmp_path, bad_batch, message):
    with pytest.raises(ValueError, match=message):
        stream_embeddings_to_h5(
            [("p1", "MA", "g1")],
            tmp_path / "bad.h5",
            embedding_dim=2,
            model_name="test-model",
            batch_size=1,
            embed_batch=lambda _sequences: bad_batch,
        )
