"""
local_constraints.py

Geometric local constraints for NanoAmorphizer.

This module provides species-sensitive checks used to reject atomic
displacements that generate severe local overlaps or excessive
distortion of the parent first-neighbor network.

The constraints are geometric sanity checks. They are not an
energy-minimization method and do not define equilibrium bond lengths.
"""

from __future__ import annotations

import numpy as np

from ase import Atoms
from ase.data import atomic_numbers, covalent_radii


# Provisional geometric safety factors.
# These values are not physical constants or material-independent
# energetic criteria. They are documented to preserve reproducibility
# and may be refined in future validated versions.
DEFAULT_MIN_DISTANCE_SCALE = 0.70
DEFAULT_PARENT_MIN_SCALE = 0.85
DEFAULT_PARENT_MAX_SCALE = 1.20
DEFAULT_NEIGHBOR_CUTOFF_SCALE = 1.20
DEFAULT_MAX_ANGLE_DEVIATION = 20.0
DEFAULT_MAX_ATTEMPTS = 50


def _validate_distance_scale(
    min_distance_scale: float,
) -> None:
    """
    Validate the anti-overlap scaling factor.
    """

    if min_distance_scale <= 0:
        raise ValueError(
            "min_distance_scale must be greater than 0."
        )

    if min_distance_scale > 1.0:
        raise ValueError(
            "min_distance_scale must not exceed 1.0 "
            "for the hard-overlap filter."
        )


def validate_parent_distance_scales(
    parent_min_scale: float,
    parent_max_scale: float,
) -> None:
    """
    Validate the allowed distortion interval for parent first neighbors.
    """

    if parent_min_scale <= 0:
        raise ValueError(
            "parent_min_scale must be greater than 0."
        )

    if parent_max_scale <= 0:
        raise ValueError(
            "parent_max_scale must be greater than 0."
        )

    if parent_min_scale >= 1.0:
        raise ValueError(
            "parent_min_scale must be smaller than 1.0."
        )

    if parent_max_scale <= 1.0:
        raise ValueError(
            "parent_max_scale must be greater than 1.0."
        )

    if parent_min_scale >= parent_max_scale:
        raise ValueError(
            "parent_min_scale must be smaller than "
            "parent_max_scale."
        )


def validate_max_attempts(
    max_attempts: int,
) -> None:
    """
    Validate the maximum number of displacement attempts.
    """

    if not isinstance(max_attempts, int):
        raise TypeError(
            "max_attempts must be an integer."
        )

    if max_attempts <= 0:
        raise ValueError(
            "max_attempts must be greater than 0."
        )


def _get_covalent_radius_from_atomic_number(
    atomic_number: int,
) -> float:
    """
    Return the ASE covalent radius for one element.
    """

    radius = float(covalent_radii[atomic_number])

    if not np.isfinite(radius) or radius <= 0:
        raise ValueError(
            f"No valid covalent radius is available "
            f"for atomic number {atomic_number}."
        )

    return radius


def get_minimum_pair_distance(
    element_a: str,
    element_b: str,
    min_distance_scale: float = DEFAULT_MIN_DISTANCE_SCALE,
) -> float:
    """
    Estimate a hard minimum separation for two chemical species.

    The limit is calculated from ASE covalent radii:

        d_min = scale * (r_cov_A + r_cov_B)

    This value is used only as an anti-overlap threshold.
    """

    _validate_distance_scale(min_distance_scale)

    try:
        atomic_number_a = atomic_numbers[element_a]
        atomic_number_b = atomic_numbers[element_b]

    except KeyError as exc:
        raise ValueError(
            f"Unknown chemical symbol: {exc.args[0]}"
        ) from exc

    radius_a = _get_covalent_radius_from_atomic_number(
        atomic_number_a
    )

    radius_b = _get_covalent_radius_from_atomic_number(
        atomic_number_b
    )

    return min_distance_scale * (
        radius_a + radius_b
    )


def get_overlap_violations(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    min_distance_scale: float = DEFAULT_MIN_DISTANCE_SCALE,
) -> dict:
    """
    Check a proposed atomic position against all other atoms.

    The proposed coordinate is compared with every other atom using
    species-sensitive minimum distances derived from ASE covalent radii.
    """

    _validate_distance_scale(min_distance_scale)

    total_atoms = len(atoms)

    if total_atoms == 0:
        raise ValueError(
            "Cannot check overlaps in an empty structure."
        )

    if atom_index < 0 or atom_index >= total_atoms:
        raise IndexError(
            "atom_index is outside the valid atom range."
        )

    proposed_position = np.asarray(
        proposed_position,
        dtype=float,
    )

    if proposed_position.shape != (3,):
        raise ValueError(
            "proposed_position must contain exactly "
            "three Cartesian coordinates."
        )

    positions = atoms.get_positions()
    atomic_numbers_array = atoms.get_atomic_numbers()

    radii = covalent_radii[
        atomic_numbers_array
    ].astype(float)

    if np.any(~np.isfinite(radii)) or np.any(radii <= 0):
        raise ValueError(
            "At least one element does not have a valid "
            "ASE covalent radius."
        )

    candidate_radius = radii[atom_index]

    displacement_vectors = (
        positions - proposed_position
    )

    distances = np.linalg.norm(
        displacement_vectors,
        axis=1,
    )

    minimum_allowed = (
        min_distance_scale
        * (candidate_radius + radii)
    )

    valid_other_atom = np.ones(
        total_atoms,
        dtype=bool,
    )

    valid_other_atom[atom_index] = False

    violation_mask = (
        valid_other_atom
        & (distances < minimum_allowed)
    )

    violating_indices = np.flatnonzero(
        violation_mask
    )

    symbols = atoms.get_chemical_symbols()

    return {
        "has_overlap": bool(
            len(violating_indices) > 0
        ),
        "violating_indices": violating_indices,
        "distances_A": distances[
            violating_indices
        ],
        "minimum_distances_A": minimum_allowed[
            violating_indices
        ],
        "candidate_symbol": symbols[
            atom_index
        ],
        "other_symbols": [
            symbols[index]
            for index in violating_indices
        ],
    }


def check_overlap(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    min_distance_scale: float = DEFAULT_MIN_DISTANCE_SCALE,
) -> bool:
    """
    Return True when a proposed position passes the anti-overlap check.
    """

    result = get_overlap_violations(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        min_distance_scale=min_distance_scale,
    )

    return not result["has_overlap"]


def build_parent_neighbor_reference(
    neighbor_network,
    n_atoms: int,
) -> list:
    """
    Convert the parent first-neighbor network into an atom-centered map.

    Parameters
    ----------
    neighbor_network : list
        Network generated by neighbors.build_neighbor_network().
        Each entry must have the form:

            (atom_i, atom_j, reference_distance)

    n_atoms : int
        Number of atoms in the nanoparticle.

    Returns
    -------
    list
        For each atom, a list of:

            (neighbor_index, reference_distance)
    """

    if n_atoms < 0:
        raise ValueError(
            "n_atoms must not be negative."
        )

    reference = [
        [] for _ in range(n_atoms)
    ]

    for i, j, reference_distance in neighbor_network:

        i = int(i)
        j = int(j)
        reference_distance = float(reference_distance)

        if i < 0 or i >= n_atoms or j < 0 or j >= n_atoms:
            raise IndexError(
                "Neighbor-network index is outside "
                "the valid atom range."
            )

        if reference_distance <= 0:
            raise ValueError(
                "Reference neighbor distances must be "
                "greater than 0 Å."
            )

        reference[i].append(
            (j, reference_distance)
        )

        reference[j].append(
            (i, reference_distance)
        )

    return reference



def get_allowed_pair_types(
    atoms: Atoms,
    neighbor_network,
) -> set:
    """
    Derive the chemical pair types present in the parent
    first-neighbor network.

    A-B and B-A are treated as the same pair type.
    """

    symbols = atoms.get_chemical_symbols()
    allowed_pair_types = set()

    for i, j, _ in neighbor_network:
        pair_type = tuple(
            sorted(
                (
                    symbols[int(i)],
                    symbols[int(j)],
                )
            )
        )
        allowed_pair_types.add(pair_type)

    return allowed_pair_types


def check_allowed_neighbor_types(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    allowed_pair_types,
    neighbor_cutoff_scale: float = DEFAULT_NEIGHBOR_CUTOFF_SCALE,
) -> bool:
    """
    Reject a proposed position if it creates a first-neighbor
    chemical pair type absent from the parent neighbor network.
    """

    if neighbor_cutoff_scale <= 0:
        raise ValueError(
            "neighbor_cutoff_scale must be greater than 0."
        )

    total_atoms = len(atoms)

    if atom_index < 0 or atom_index >= total_atoms:
        raise IndexError(
            "atom_index is outside the valid atom range."
        )

    proposed_position = np.asarray(
        proposed_position,
        dtype=float,
    )

    if proposed_position.shape != (3,):
        raise ValueError(
            "proposed_position must contain exactly "
            "three Cartesian coordinates."
        )

    positions = atoms.get_positions()
    atomic_numbers_array = atoms.get_atomic_numbers()
    symbols = atoms.get_chemical_symbols()

    radii = covalent_radii[
        atomic_numbers_array
    ].astype(float)

    if np.any(~np.isfinite(radii)) or np.any(radii <= 0):
        raise ValueError(
            "At least one element does not have a valid "
            "ASE covalent radius."
        )

    candidate_radius = radii[atom_index]
    candidate_symbol = symbols[atom_index]

    distances = np.linalg.norm(
        positions - proposed_position,
        axis=1,
    )

    neighbor_limits = (
        neighbor_cutoff_scale
        * (candidate_radius + radii)
    )

    for other_index in range(total_atoms):

        if other_index == atom_index:
            continue

        if distances[other_index] > neighbor_limits[other_index]:
            continue

        pair_type = tuple(
            sorted(
                (
                    candidate_symbol,
                    symbols[other_index],
                )
            )
        )

        if pair_type not in allowed_pair_types:
            return False

    return True


def validate_max_angle_deviation(
    max_angle_deviation: float,
) -> None:
    """
    Validate the maximum allowed angular deviation in degrees.
    """

    if max_angle_deviation <= 0:
        raise ValueError(
            "max_angle_deviation must be greater than 0 degrees."
        )

    if max_angle_deviation >= 180:
        raise ValueError(
            "max_angle_deviation must be smaller than 180 degrees."
        )


def _calculate_angle_degrees(
    position_j,
    position_i,
    position_k,
) -> float:
    """
    Calculate angle j-i-k in degrees, with atom i as the vertex.
    """

    vector_ij = np.asarray(position_j, dtype=float) - np.asarray(
        position_i,
        dtype=float,
    )
    vector_ik = np.asarray(position_k, dtype=float) - np.asarray(
        position_i,
        dtype=float,
    )

    norm_ij = np.linalg.norm(vector_ij)
    norm_ik = np.linalg.norm(vector_ik)

    if norm_ij <= 0 or norm_ik <= 0:
        raise ValueError(
            "Cannot calculate an angle from a zero-length vector."
        )

    cosine = np.dot(vector_ij, vector_ik) / (
        norm_ij * norm_ik
    )

    cosine = np.clip(cosine, -1.0, 1.0)

    return float(
        np.degrees(
            np.arccos(cosine)
        )
    )


def build_parent_angle_reference(
    atoms: Atoms,
    parent_neighbor_reference,
) -> list:
    """
    Build the angular reference from the parent crystalline structure.

    For each central atom i and every unique pair of its parent
    first neighbors j and k, the reference angle j-i-k is stored.

    The returned structure is indexed by atom. Each atom therefore
    has access to every parent angle in which it participates,
    whether as the central atom or as one of the two outer atoms.
    """

    total_atoms = len(atoms)

    if len(parent_neighbor_reference) != total_atoms:
        raise ValueError(
            "parent_neighbor_reference must contain one "
            "entry for every atom."
        )

    positions = atoms.get_positions()

    angle_reference = [
        [] for _ in range(total_atoms)
    ]

    for center_index in range(total_atoms):

        neighbor_indices = [
            int(neighbor_index)
            for neighbor_index, _ in (
                parent_neighbor_reference[center_index]
            )
        ]

        for first in range(
            len(neighbor_indices) - 1
        ):
            for second in range(
                first + 1,
                len(neighbor_indices),
            ):

                atom_j = neighbor_indices[first]
                atom_k = neighbor_indices[second]

                reference_angle = _calculate_angle_degrees(
                    positions[atom_j],
                    positions[center_index],
                    positions[atom_k],
                )

                angle_record = (
                    atom_j,
                    center_index,
                    atom_k,
                    reference_angle,
                )

                angle_reference[atom_j].append(
                    angle_record
                )
                angle_reference[center_index].append(
                    angle_record
                )
                angle_reference[atom_k].append(
                    angle_record
                )

    return angle_reference


def get_parent_angle_violations(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    parent_angle_reference,
    max_angle_deviation: float = DEFAULT_MAX_ANGLE_DEVIATION,
) -> dict:
    """
    Check all parent angles involving a proposed atomic position.

    A parent angle is accepted when:

        abs(theta_new - theta_parent) <= max_angle_deviation
    """

    validate_max_angle_deviation(
        max_angle_deviation
    )

    total_atoms = len(atoms)

    if atom_index < 0 or atom_index >= total_atoms:
        raise IndexError(
            "atom_index is outside the valid atom range."
        )

    if len(parent_angle_reference) != total_atoms:
        raise ValueError(
            "parent_angle_reference must contain one "
            "entry for every atom."
        )

    proposed_position = np.asarray(
        proposed_position,
        dtype=float,
    )

    if proposed_position.shape != (3,):
        raise ValueError(
            "proposed_position must contain exactly "
            "three Cartesian coordinates."
        )

    positions = atoms.get_positions()

    violating_angles = []

    for (
        atom_j,
        center_index,
        atom_k,
        reference_angle,
    ) in parent_angle_reference[atom_index]:

        position_j = (
            proposed_position
            if atom_j == atom_index
            else positions[atom_j]
        )

        position_i = (
            proposed_position
            if center_index == atom_index
            else positions[center_index]
        )

        position_k = (
            proposed_position
            if atom_k == atom_index
            else positions[atom_k]
        )

        proposed_angle = _calculate_angle_degrees(
            position_j,
            position_i,
            position_k,
        )

        deviation = abs(
            proposed_angle - reference_angle
        )

        if deviation > max_angle_deviation:
            violating_angles.append(
                {
                    "atom_j": int(atom_j),
                    "center_atom": int(center_index),
                    "atom_k": int(atom_k),
                    "reference_angle_deg": float(reference_angle),
                    "proposed_angle_deg": float(proposed_angle),
                    "deviation_deg": float(deviation),
                }
            )

    return {
        "has_violation": bool(
            len(violating_angles) > 0
        ),
        "violations": violating_angles,
    }


def check_parent_angles(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    parent_angle_reference,
    max_angle_deviation: float = DEFAULT_MAX_ANGLE_DEVIATION,
) -> bool:
    """
    Return True when all parent angles involving the proposed atom
    remain within the allowed angular deviation.
    """

    result = get_parent_angle_violations(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        parent_angle_reference=parent_angle_reference,
        max_angle_deviation=max_angle_deviation,
    )

    return not result["has_violation"]

def get_parent_neighbor_violations(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    parent_neighbor_reference,
    parent_min_scale: float = DEFAULT_PARENT_MIN_SCALE,
    parent_max_scale: float = DEFAULT_PARENT_MAX_SCALE,
) -> dict:
    """
    Check distortion of the parent first-neighbor network.

    For each parent first neighbor of atom_index, the proposed distance
    must remain inside:

        parent_min_scale * d0 <= d_new <= parent_max_scale * d0

    where d0 is the corresponding distance in the parent crystalline
    nanoparticle.
    """

    validate_parent_distance_scales(
        parent_min_scale,
        parent_max_scale,
    )

    total_atoms = len(atoms)

    if atom_index < 0 or atom_index >= total_atoms:
        raise IndexError(
            "atom_index is outside the valid atom range."
        )

    if len(parent_neighbor_reference) != total_atoms:
        raise ValueError(
            "parent_neighbor_reference must contain one "
            "entry for every atom."
        )

    proposed_position = np.asarray(
        proposed_position,
        dtype=float,
    )

    if proposed_position.shape != (3,):
        raise ValueError(
            "proposed_position must contain exactly "
            "three Cartesian coordinates."
        )

    positions = atoms.get_positions()

    violating_indices = []
    reference_distances = []
    proposed_distances = []
    minimum_distances = []
    maximum_distances = []

    for neighbor_index, reference_distance in (
        parent_neighbor_reference[atom_index]
    ):

        neighbor_position = positions[
            neighbor_index
        ]

        proposed_distance = float(
            np.linalg.norm(
                proposed_position - neighbor_position
            )
        )

        minimum_distance = (
            parent_min_scale
            * reference_distance
        )

        maximum_distance = (
            parent_max_scale
            * reference_distance
        )

        if (
            proposed_distance < minimum_distance
            or proposed_distance > maximum_distance
        ):
            violating_indices.append(
                int(neighbor_index)
            )
            reference_distances.append(
                float(reference_distance)
            )
            proposed_distances.append(
                proposed_distance
            )
            minimum_distances.append(
                minimum_distance
            )
            maximum_distances.append(
                maximum_distance
            )

    symbols = atoms.get_chemical_symbols()

    return {
        "has_violation": bool(
            len(violating_indices) > 0
        ),
        "violating_indices": np.asarray(
            violating_indices,
            dtype=int,
        ),
        "reference_distances_A": np.asarray(
            reference_distances,
            dtype=float,
        ),
        "proposed_distances_A": np.asarray(
            proposed_distances,
            dtype=float,
        ),
        "minimum_distances_A": np.asarray(
            minimum_distances,
            dtype=float,
        ),
        "maximum_distances_A": np.asarray(
            maximum_distances,
            dtype=float,
        ),
        "candidate_symbol": symbols[
            atom_index
        ],
        "other_symbols": [
            symbols[index]
            for index in violating_indices
        ],
    }


def check_parent_neighbors(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    parent_neighbor_reference,
    parent_min_scale: float = DEFAULT_PARENT_MIN_SCALE,
    parent_max_scale: float = DEFAULT_PARENT_MAX_SCALE,
) -> bool:
    """
    Return True when all parent first-neighbor distances remain
    inside the allowed distortion interval.
    """

    result = get_parent_neighbor_violations(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        parent_neighbor_reference=parent_neighbor_reference,
        parent_min_scale=parent_min_scale,
        parent_max_scale=parent_max_scale,
    )

    return not result["has_violation"]


def evaluate_local_constraints(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    parent_neighbor_reference,
    parent_angle_reference,
    allowed_pair_types,
    min_distance_scale: float = DEFAULT_MIN_DISTANCE_SCALE,
    parent_min_scale: float = DEFAULT_PARENT_MIN_SCALE,
    parent_max_scale: float = DEFAULT_PARENT_MAX_SCALE,
    neighbor_cutoff_scale: float = DEFAULT_NEIGHBOR_CUTOFF_SCALE,
    max_angle_deviation: float = DEFAULT_MAX_ANGLE_DEVIATION,
) -> dict:
    """
    Evaluate all active local geometric constraints.

    Returns
    -------
    dict
        {
            "accepted": bool,
            "reason": str or None
        }

    Notes
    -----
    The rejection reason is the first failed constraint in the
    evaluation sequence. It is therefore a diagnostic category,
    not an exhaustive list of every condition the proposal might
    violate simultaneously.
    """

    if not check_overlap(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        min_distance_scale=min_distance_scale,
    ):
        return {
            "accepted": False,
            "reason": "overlap",
        }

    if not check_parent_neighbors(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        parent_neighbor_reference=parent_neighbor_reference,
        parent_min_scale=parent_min_scale,
        parent_max_scale=parent_max_scale,
    ):
        return {
            "accepted": False,
            "reason": "parent_distance",
        }

    if not check_allowed_neighbor_types(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        allowed_pair_types=allowed_pair_types,
        neighbor_cutoff_scale=neighbor_cutoff_scale,
    ):
        return {
            "accepted": False,
            "reason": "pair_type",
        }

    if not check_parent_angles(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        parent_angle_reference=parent_angle_reference,
        max_angle_deviation=max_angle_deviation,
    ):
        return {
            "accepted": False,
            "reason": "angle",
        }

    return {
        "accepted": True,
        "reason": None,
    }


def check_local_constraints(
    atoms: Atoms,
    atom_index: int,
    proposed_position,
    parent_neighbor_reference,
    parent_angle_reference,
    allowed_pair_types,
    min_distance_scale: float = DEFAULT_MIN_DISTANCE_SCALE,
    parent_min_scale: float = DEFAULT_PARENT_MIN_SCALE,
    parent_max_scale: float = DEFAULT_PARENT_MAX_SCALE,
    neighbor_cutoff_scale: float = DEFAULT_NEIGHBOR_CUTOFF_SCALE,
    max_angle_deviation: float = DEFAULT_MAX_ANGLE_DEVIATION,
) -> bool:
    """
    Return True when a proposed position passes every active
    local geometric constraint.
    """

    result = evaluate_local_constraints(
        atoms=atoms,
        atom_index=atom_index,
        proposed_position=proposed_position,
        parent_neighbor_reference=parent_neighbor_reference,
        parent_angle_reference=parent_angle_reference,
        allowed_pair_types=allowed_pair_types,
        min_distance_scale=min_distance_scale,
        parent_min_scale=parent_min_scale,
        parent_max_scale=parent_max_scale,
        neighbor_cutoff_scale=neighbor_cutoff_scale,
        max_angle_deviation=max_angle_deviation,
    )

    return result["accepted"]
