"""
geometry.py

Geometric construction and cleanup utilities for NanoAmorphizer.
"""

import numpy as np

from neighbors import (
    build_neighbor_network,
    get_coordination_numbers,
)


def get_particle_center(atoms):
    """
    Return the geometric center of the atomic coordinates.
    """

    if len(atoms) == 0:
        raise ValueError(
            "Cannot determine the center of an empty structure."
        )

    positions = atoms.get_positions()

    return positions.mean(axis=0)


def crop_sphere(
    atoms,
    radius,
    center=None,
):
    """
    Crop a finite spherical nanoparticle from an atomic structure.

    Returns
    -------
    atoms_sphere : ase.Atoms
        Non-periodic spherical particle.

    center : numpy.ndarray
        Cartesian center used for the spherical crop.
    """

    if radius <= 0:
        raise ValueError(
            "Particle radius must be greater than 0."
        )

    positions = atoms.get_positions()

    if center is None:
        center = get_particle_center(atoms)
    else:
        center = np.asarray(
            center,
            dtype=float,
        )

        if center.shape != (3,):
            raise ValueError(
                "center must contain exactly three "
                "Cartesian coordinates."
            )

    distances = np.linalg.norm(
        positions - center,
        axis=1,
    )

    mask = distances <= radius

    atoms_sphere = atoms[mask].copy()

    # NanoAmorphizer treats nanoparticles as finite,
    # non-periodic structures.
    atoms_sphere.set_pbc(False)

    return atoms_sphere, center.copy()


def remove_isolated_atoms(
    atoms,
    cutoff_scale=1.20,
):
    """
    Remove atoms with zero first-neighbor coordination.

    Isolation is evaluated with the same first-neighbor definition
    used by neighbors.build_neighbor_network(), based on ASE covalent
    radii and the selected cutoff scale.

    Only atoms with coordination number equal to zero are removed.
    Under-coordinated surface atoms with one or more neighbors are
    retained.

    Returns
    -------
    cleaned_atoms : ase.Atoms
        Copy of the structure without isolated atoms.

    removed_count : int
        Number of atoms removed.
    """

    if len(atoms) == 0:
        return atoms.copy(), 0

    if cutoff_scale <= 0:
        raise ValueError(
            "cutoff_scale must be greater than 0."
        )

    neighbor_network = build_neighbor_network(
        atoms,
        cutoff_scale=cutoff_scale,
    )

    coordination = get_coordination_numbers(
        atoms,
        neighbor_network,
    )

    isolated_mask = coordination == 0
    removed_count = int(
        np.count_nonzero(isolated_mask)
    )

    if removed_count == 0:
        return atoms.copy(), 0

    if removed_count == len(atoms):
        raise ValueError(
            "All atoms were classified as isolated. "
            "Check the neighbor cutoff or input structure."
        )

    keep_mask = ~isolated_mask

    cleaned_atoms = atoms[
        keep_mask
    ].copy()

    cleaned_atoms.set_pbc(False)

    return cleaned_atoms, removed_count
