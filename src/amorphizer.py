import numpy as np

from shell import get_shell_mask, get_shell_weights
from neighbors import build_neighbor_network
from local_constraints import (
    DEFAULT_MAX_ATTEMPTS,
    DEFAULT_MIN_DISTANCE_SCALE,
    DEFAULT_PARENT_MIN_SCALE,
    DEFAULT_PARENT_MAX_SCALE,
    DEFAULT_MAX_ANGLE_DEVIATION,
    build_parent_neighbor_reference,
    build_parent_angle_reference,
    get_allowed_pair_types,
    evaluate_local_constraints,
    validate_max_attempts,
    validate_parent_distance_scales,
    validate_max_angle_deviation,
)


def amorphize_shell(
    atoms,
    radius,
    shell_thickness,
    disorder_strength,
    seed=None,
    center=None,
    min_distance_scale=DEFAULT_MIN_DISTANCE_SCALE,
    parent_min_scale=DEFAULT_PARENT_MIN_SCALE,
    parent_max_scale=DEFAULT_PARENT_MAX_SCALE,
    max_angle_deviation=DEFAULT_MAX_ANGLE_DEVIATION,
    max_attempts=DEFAULT_MAX_ATTEMPTS,
    return_report=False,
):
    """
    Apply progressively increasing stochastic displacement to the
    surface shell of a spherical nanoparticle.

    The crystalline core is preserved. Atomic displacement increases
    from the core-shell boundary toward the nominal particle surface
    according to the radial weights defined in shell.py.

    Every proposed shell displacement is checked against local
    geometric constraints derived from the parent structure:

    1. a species-sensitive hard anti-overlap threshold,
    2. the allowed distortion interval of the parent first-neighbor network,
    3. allowed first-neighbor chemical pair types, and
    4. allowed deviation from parent first-neighbor angles.

    Invalid proposals are rejected and regenerated up to max_attempts
    times. If no valid proposal is found, the atom remains at its
    original position.

    Parameters
    ----------
    atoms : ase.Atoms
        Finite crystalline nanoparticle.

    radius : float
        Nominal particle radius in angstroms.

    shell_thickness : float
        Thickness of the surface region affected by disorder,
        in angstroms.

    disorder_strength : float
        Maximum three-dimensional RMS displacement magnitude at the
        nominal particle surface, in angstroms.

        For atoms at the particle surface:

            RMS(|dr|) = disorder_strength

        For atoms inside the shell:

            RMS(|dr|) = disorder_strength * radial_weight

    seed : int or None, optional
        Seed used by the NumPy random number generator. Using the same
        seed and identical input parameters produces the same stochastic
        displacement sequence.

    center : array-like, optional
        Cartesian coordinates of the particle center. The same center
        used during spherical cropping should be supplied to maintain
        consistent shell geometry.

    min_distance_scale : float, optional
        Species-sensitive hard-overlap scale used by local_constraints.py.
        The threshold is based on the sum of ASE covalent radii.
        This is a geometric safety parameter, not an equilibrium
        bond-length model.

    parent_min_scale : float, optional
        Minimum allowed fraction of each parent first-neighbor distance.

    parent_max_scale : float, optional
        Maximum allowed fraction of each parent first-neighbor distance.

    max_angle_deviation : float, optional
        Maximum allowed deviation, in degrees, from each parent
        first-neighbor angle involving a displaced atom.

    max_attempts : int, optional
        Maximum number of displacement proposals tested for each shell
        atom before retaining its original position.

    return_report : bool, optional
        If True, return both the amorphized structure and a diagnostic
        report describing proposal acceptance and rejection statistics.

    Returns
    -------
    ase.Atoms or tuple
        By default, returns the amorphized nanoparticle. If
        return_report=True, returns:

            (amorphized_atoms, report)

    Notes
    -----
    This procedure generates geometrically disordered structures.
    It does not perform energy minimization or simulate a physical
    amorphization pathway.
    """

    amorphized = atoms.copy()

    report = {
        "selected_shell_atoms": 0,
        "zero_weight_atoms": 0,
        "atoms_attempted": 0,
        "atoms_displaced": 0,
        "atoms_retained_after_max_attempts": 0,
        "total_proposals": 0,
        "accepted_proposals": 0,
        "rejected_proposals": 0,
        "rejected_overlap": 0,
        "rejected_parent_distance": 0,
        "rejected_pair_type": 0,
        "rejected_angle": 0,
        "proposal_acceptance_rate_percent": 0.0,
        "shell_displacement_rate_percent": 0.0,
        "mean_proposals_per_attempted_atom": 0.0,
        "max_attempts_per_atom": int(max_attempts),
    }

    if shell_thickness <= 0 or disorder_strength <= 0:
        if return_report:
            return amorphized, report
        return amorphized

    validate_max_attempts(max_attempts)

    validate_parent_distance_scales(
        parent_min_scale,
        parent_max_scale,
    )

    validate_max_angle_deviation(
        max_angle_deviation
    )

    # Build the first-neighbor reference from the original crystalline
    # nanoparticle before any atomic displacement is applied.
    parent_neighbor_network = build_neighbor_network(
        atoms
    )

    parent_neighbor_reference = build_parent_neighbor_reference(
        parent_neighbor_network,
        n_atoms=len(atoms),
    )

    parent_angle_reference = build_parent_angle_reference(
        atoms,
        parent_neighbor_reference,
    )

    allowed_pair_types = get_allowed_pair_types(
        atoms,
        parent_neighbor_network,
    )

    shell_mask = get_shell_mask(
        atoms,
        radius=radius,
        shell_thickness=shell_thickness,
        center=center,
    )

    weights = get_shell_weights(
        atoms,
        radius=radius,
        shell_thickness=shell_thickness,
        center=center,
    )

    shell_indices = np.flatnonzero(shell_mask)
    report["selected_shell_atoms"] = int(len(shell_indices))

    rng = np.random.default_rng(seed)

    original_positions = atoms.get_positions()

    surface_sigma_xyz = disorder_strength / np.sqrt(3.0)

    for atom_index in shell_indices:

        radial_weight = float(weights[atom_index])

        if radial_weight <= 0.0:
            report["zero_weight_atoms"] += 1
            continue

        report["atoms_attempted"] += 1

        original_position = original_positions[
            atom_index
        ].copy()

        accepted = False

        for _ in range(max_attempts):

            report["total_proposals"] += 1

            displacement = rng.normal(
                loc=0.0,
                scale=surface_sigma_xyz,
                size=3,
            )

            displacement *= radial_weight

            proposed_position = (
                original_position + displacement
            )

            constraint_result = (
                evaluate_local_constraints(
                    atoms=amorphized,
                    atom_index=int(atom_index),
                    proposed_position=proposed_position,
                    parent_neighbor_reference=parent_neighbor_reference,
                    parent_angle_reference=parent_angle_reference,
                    allowed_pair_types=allowed_pair_types,
                    min_distance_scale=min_distance_scale,
                    parent_min_scale=parent_min_scale,
                    parent_max_scale=parent_max_scale,
                    max_angle_deviation=max_angle_deviation,
                )
            )

            if constraint_result["accepted"]:

                amorphized.positions[
                    atom_index
                ] = proposed_position

                report["accepted_proposals"] += 1
                report["atoms_displaced"] += 1
                accepted = True
                break

            report["rejected_proposals"] += 1

            reason = constraint_result["reason"]

            if reason == "overlap":
                report["rejected_overlap"] += 1

            elif reason == "parent_distance":
                report["rejected_parent_distance"] += 1

            elif reason == "pair_type":
                report["rejected_pair_type"] += 1

            elif reason == "angle":
                report["rejected_angle"] += 1

        if not accepted:
            report[
                "atoms_retained_after_max_attempts"
            ] += 1

    if report["total_proposals"] > 0:
        report[
            "proposal_acceptance_rate_percent"
        ] = (
            100.0
            * report["accepted_proposals"]
            / report["total_proposals"]
        )

    if report["selected_shell_atoms"] > 0:
        report[
            "shell_displacement_rate_percent"
        ] = (
            100.0
            * report["atoms_displaced"]
            / report["selected_shell_atoms"]
        )

    if report["atoms_attempted"] > 0:
        report[
            "mean_proposals_per_attempted_atom"
        ] = (
            report["total_proposals"]
            / report["atoms_attempted"]
        )

    if return_report:
        return amorphized, report

    return amorphized

