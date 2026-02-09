"""
Transmembrane topology prediction using pyTMHMM.

This module provides integration with pyTMHMM for predicting transmembrane
helices and deriving topology-based predicates.

Requires the 'topology' extra: pip install sharur[topology]
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import re


# Try to import pyTMHMM, but don't fail if not installed
_PYTMHMM_AVAILABLE = False
try:
    import pyTMHMM
    _PYTMHMM_AVAILABLE = True
except ImportError:
    pyTMHMM = None  # type: ignore


def is_available() -> bool:
    """Check if pyTMHMM is available."""
    return _PYTMHMM_AVAILABLE


@dataclass
class TopologyPrediction:
    """Result of transmembrane topology prediction."""

    sequence_length: int
    num_tm_helices: int
    tm_helix_positions: list[tuple[int, int]]  # (start, end) positions
    topology: str  # Full topology string (i/o/M annotation per residue)
    n_terminus_location: str  # 'inside' or 'outside'

    @property
    def is_transmembrane(self) -> bool:
        """Whether the protein has any predicted TM helices."""
        return self.num_tm_helices > 0

    @property
    def is_single_pass(self) -> bool:
        """Whether the protein has exactly one TM helix."""
        return self.num_tm_helices == 1

    @property
    def is_multi_pass(self) -> bool:
        """Whether the protein has 2+ TM helices."""
        return self.num_tm_helices >= 2

    @property
    def is_polytopic(self) -> bool:
        """Whether the protein has 4+ TM helices (polytopic membrane protein)."""
        return self.num_tm_helices >= 4


def predict_topology(sequence: str, compute_posterior: bool = False) -> Optional[TopologyPrediction]:
    """
    Predict transmembrane topology for a protein sequence.

    Args:
        sequence: Amino acid sequence (single-letter code)
        compute_posterior: Whether to compute posterior probabilities (slower)

    Returns:
        TopologyPrediction object, or None if pyTMHMM is not available

    Raises:
        ValueError: If sequence is empty or contains invalid characters
    """
    if not _PYTMHMM_AVAILABLE:
        return None

    if not sequence:
        raise ValueError("Sequence cannot be empty")

    # Clean sequence - remove whitespace and convert to uppercase
    sequence = re.sub(r'\s', '', sequence.upper())

    # Validate sequence contains only valid amino acids
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    invalid = set(sequence) - valid_aa
    if invalid:
        # Replace invalid characters with X (unknown)
        sequence = ''.join(c if c in valid_aa else 'X' for c in sequence)

    # Run pyTMHMM prediction
    annotation = pyTMHMM.predict(sequence, compute_posterior=compute_posterior)

    # If posterior was requested, annotation is a tuple
    if compute_posterior:
        annotation = annotation[0]

    # Parse the annotation string to extract TM helix positions
    # Annotation uses: 'M' = transmembrane, 'i' = inside, 'o' = outside
    tm_positions = _parse_tm_positions(annotation)

    # Determine N-terminus location
    n_term = 'inside' if annotation[0] == 'i' else 'outside' if annotation[0] == 'o' else 'membrane'

    return TopologyPrediction(
        sequence_length=len(sequence),
        num_tm_helices=len(tm_positions),
        tm_helix_positions=tm_positions,
        topology=annotation,
        n_terminus_location=n_term,
    )


def _parse_tm_positions(annotation: str) -> list[tuple[int, int]]:
    """
    Parse TM helix positions from annotation string.

    Args:
        annotation: String of 'M', 'i', 'o' characters (one per residue)

    Returns:
        List of (start, end) tuples for each TM helix (1-indexed)
    """
    positions = []
    in_helix = False
    start = 0

    for i, char in enumerate(annotation):
        if char == 'M' and not in_helix:
            # Start of new TM helix
            in_helix = True
            start = i + 1  # 1-indexed
        elif char != 'M' and in_helix:
            # End of TM helix
            in_helix = False
            positions.append((start, i))  # end is exclusive, so i (not i+1)

    # Handle case where sequence ends in TM helix
    if in_helix:
        positions.append((start, len(annotation)))

    return positions


def get_topology_predicates(prediction: TopologyPrediction) -> list[str]:
    """
    Get predicate IDs from a topology prediction.

    Args:
        prediction: TopologyPrediction object

    Returns:
        List of predicate IDs
    """
    predicates = []

    if prediction.is_transmembrane:
        predicates.append("transmembrane_predicted")

        if prediction.is_single_pass:
            predicates.append("single_pass_membrane")
        elif prediction.is_multi_pass:
            predicates.append("multi_pass_membrane")

            if prediction.is_polytopic:
                predicates.append("polytopic_membrane")

        # Topology orientation
        if prediction.n_terminus_location == 'inside':
            predicates.append("n_in_topology")
        elif prediction.n_terminus_location == 'outside':
            predicates.append("n_out_topology")
    else:
        predicates.append("soluble_predicted")

    return predicates


def predict_topology_batch(
    sequences: dict[str, str],
) -> dict[str, TopologyPrediction]:
    """
    Predict topology for multiple sequences.

    Args:
        sequences: Dict mapping protein_id -> sequence

    Returns:
        Dict mapping protein_id -> TopologyPrediction
        (only proteins that could be predicted are included)
    """
    if not _PYTMHMM_AVAILABLE:
        return {}

    results = {}
    for protein_id, sequence in sequences.items():
        try:
            pred = predict_topology(sequence)
            if pred is not None:
                results[protein_id] = pred
        except (ValueError, Exception):
            # Skip proteins with invalid sequences
            continue

    return results


__all__ = [
    "is_available",
    "TopologyPrediction",
    "predict_topology",
    "get_topology_predicates",
    "predict_topology_batch",
]
