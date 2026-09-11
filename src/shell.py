"""
shell.py

Geometric definition of the core and surface shell in spherical nanoparticles.

This module identifies atoms belonging to the crystalline core, the selected
surface shell, and (when relevant) atoms located outside the nominal particle
radius.

NanoAmorphizer
"""

from __future__ import annotations

import numpy as np
from ase import Atoms


_GEOMETRIC_TOLERANCE = 1e-8


def linear_disorder_profile(
    normalized_radius: np.ndarray,
) -> np.ndarray:
    """
    Return a linear radial disorder profile.

    The input is the normalized shell coordinate:

        u = (r - (R - t)) / t

    where:
        u = 0 at the core-shell boundary
        u = 1 at the nominal particle surface

    The linear profile is:

        f(u) = u

    Values are clipped to the interval [0, 1].

    Parameters
    ----------
    normalized_radius : array-like
        Normalized radial coordinate within the shell.

    Returns
    -------
    numpy.ndarray
        Linear disorder weights in the interval [0, 1].
    """

    normalized_radius = np.asarray(
        normalized_radius,
        dtype=float,
    )

    return np.clip(
        normalized_radius,
        0.0,
        1.0,
    )


def _validate_shell_parameters(
    radius: float,
    shell_thickness: float,
) -> None:
    """
    Validate spherical shell parameters.

    Parameters
    ----------
    radius : float
        Particle radius in angstroms.

    shell_thickness : float
        Shell thickness in angstroms.

    Raises
    ------
    ValueError
        If radius or shell thickness are outside the allowed range.
    """

    if radius <= 0:
        raise ValueError("Particle radius must be greater than 0 Å.")

    if shell_thickness <= 0:
        raise ValueError("Shell thickness must be greater than 0 Å.")

    if shell_thickness >= radius:
        raise ValueError(
            "Shell thickness must be smaller than the particle radius."
        )


def get_radial_distances(
    atoms: Atoms,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Calculate the radial distance of every atom from the particle center.

    Parameters
    ----------
    atoms : ase.Atoms
        Nanoparticle atomic structure.

    center : array-like, optional
        Particle center in Cartesian coordinates. If omitted, the
        geometric center of the atomic positions is used.

    Returns
    -------
    numpy.ndarray
        Radial distance of each atom from the center, in angstroms.
    """

    if len(atoms) == 0:
        return np.array([], dtype=float)

    positions = atoms.get_positions()

    if center is None:
        center = positions.mean(axis=0)
    else:
        center = np.asarray(center, dtype=float)

        if center.shape != (3,):
            raise ValueError(
                "center must contain exactly three Cartesian coordinates."
            )

    return np.linalg.norm(positions - center, axis=1)


def get_shell_mask(
    atoms: Atoms,
    radius: float,
    shell_thickness: float,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Identify atoms belonging to the surface shell.

    The shell extends inward from the nominal particle surface:

        radius - shell_thickness <= r <= radius

    Parameters
    ----------
    atoms : ase.Atoms
        Spherical nanoparticle.

    radius : float
        Nominal particle radius in angstroms.

    shell_thickness : float
        Thickness of the surface shell in angstroms.

    center : array-like, optional
        Particle center.

    Returns
    -------
    numpy.ndarray
        Boolean mask. True values correspond to shell atoms.
    """

    _validate_shell_parameters(radius, shell_thickness)

    radial_distances = get_radial_distances(
        atoms,
        center=center,
    )

    inner_radius = radius - shell_thickness

    return (
        (radial_distances >= inner_radius - _GEOMETRIC_TOLERANCE)
        & (radial_distances <= radius + _GEOMETRIC_TOLERANCE)
    )


def get_core_mask(
    atoms: Atoms,
    radius: float,
    shell_thickness: float,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Identify atoms belonging to the crystalline core.

    The core is defined explicitly as:

        r < radius - shell_thickness

    Atoms displaced beyond the nominal particle radius are therefore
    not incorrectly reclassified as core atoms.

    Parameters
    ----------
    atoms : ase.Atoms
        Spherical nanoparticle.

    radius : float
        Nominal particle radius in angstroms.

    shell_thickness : float
        Surface shell thickness in angstroms.

    center : array-like, optional
        Particle center.

    Returns
    -------
    numpy.ndarray
        Boolean mask. True values correspond to core atoms.
    """

    _validate_shell_parameters(radius, shell_thickness)

    radial_distances = get_radial_distances(
        atoms,
        center=center,
    )

    inner_radius = radius - shell_thickness

    return radial_distances < (
        inner_radius - _GEOMETRIC_TOLERANCE
    )


def get_outside_mask(
    atoms: Atoms,
    radius: float,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Identify atoms located outside the nominal particle radius.

    This mask is mainly useful for post-amorphization validation, because
    stochastic displacements may move some shell atoms slightly beyond
    the original spherical boundary.

    Parameters
    ----------
    atoms : ase.Atoms
        Nanoparticle.

    radius : float
        Nominal particle radius in angstroms.

    center : array-like, optional
        Particle center.

    Returns
    -------
    numpy.ndarray
        Boolean mask. True values correspond to atoms with r > radius.
    """

    if radius <= 0:
        raise ValueError("Particle radius must be greater than 0 Å.")

    radial_distances = get_radial_distances(
        atoms,
        center=center,
    )

    return radial_distances > (
        radius + _GEOMETRIC_TOLERANCE
    )


def get_shell_indices(
    atoms: Atoms,
    radius: float,
    shell_thickness: float,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Return the atomic indices belonging to the surface shell.
    """

    mask = get_shell_mask(
        atoms,
        radius,
        shell_thickness,
        center=center,
    )

    return np.flatnonzero(mask)


def get_core_indices(
    atoms: Atoms,
    radius: float,
    shell_thickness: float,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Return the atomic indices belonging to the crystalline core.
    """

    mask = get_core_mask(
        atoms,
        radius,
        shell_thickness,
        center=center,
    )

    return np.flatnonzero(mask)


def get_outside_indices(
    atoms: Atoms,
    radius: float,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Return the atomic indices located outside the nominal particle radius.
    """

    mask = get_outside_mask(
        atoms,
        radius,
        center=center,
    )

    return np.flatnonzero(mask)


def get_shell_weights(
    atoms: Atoms,
    radius: float,
    shell_thickness: float,
    center: np.ndarray | None = None,
) -> np.ndarray:
    """
    Calculate normalized radial weights for atoms in the shell.

    The normalized shell coordinate is converted into a disorder weight
    through an explicit radial profile function.

    The current default profile is linear:

        f(u) = u

    where u = 0 at the core-shell boundary and u = 1 at the nominal
    particle surface.

    Core atoms receive a weight of 0. Values beyond the nominal particle
    radius are clipped to 1.

    Returns
    -------
    numpy.ndarray
        Weight for every atom in the nanoparticle.
    """

    _validate_shell_parameters(radius, shell_thickness)

    radial_distances = get_radial_distances(
        atoms,
        center=center,
    )

    inner_radius = radius - shell_thickness

    normalized_radius = (
        radial_distances - inner_radius
    ) / shell_thickness

    return linear_disorder_profile(
        normalized_radius
    )


def get_shell_info(
    atoms: Atoms,
    radius: float,
    shell_thickness: float,
    center: np.ndarray | None = None,
) -> dict:
    """
    Return a summary of the core-shell geometry.

    Core, shell and outside atoms are counted independently so that atoms
    beyond the nominal particle radius are never interpreted as core atoms.

    Returns
    -------
    dict
        Dictionary containing total, core, shell and outside atom counts,
        particle radius, shell thickness and core radius.
    """

    shell_mask = get_shell_mask(
        atoms,
        radius,
        shell_thickness,
        center=center,
    )

    core_mask = get_core_mask(
        atoms,
        radius,
        shell_thickness,
        center=center,
    )

    outside_mask = get_outside_mask(
        atoms,
        radius,
        center=center,
    )

    total_atoms = len(atoms)
    core_atoms = int(np.count_nonzero(core_mask))
    shell_atoms = int(np.count_nonzero(shell_mask))
    outside_atoms = int(np.count_nonzero(outside_mask))

    return {
        "total_atoms": total_atoms,
        "core_atoms": core_atoms,
        "shell_atoms": shell_atoms,
        "outside_atoms": outside_atoms,
        "radius_A": float(radius),
        "shell_thickness_A": float(shell_thickness),
        "core_radius_A": float(radius - shell_thickness),
    }
