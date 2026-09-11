from prompts import (
    ask_cif_path,
    ask_radius,
    ask_shell,
    ask_disorder,
    ask_seed,
)

from io_utils import load_cif
from validation import (
    validate_radius,
    validate_shell,
    validate_disorder,
)
from crystal_builder import build_supercell
from geometry import crop_sphere, remove_isolated_atoms
from export import export_xyz
from amorphizer import amorphize_shell
from shell import get_shell_info
from neighbors import (
    build_neighbor_network,
    get_neighbor_summary,
)


def get_particle_size_nm(atoms):
    positions = atoms.get_positions()
    mins = positions.min(axis=0)
    maxs = positions.max(axis=0)
    spans_nm = (maxs - mins) / 10.0

    return (
        spans_nm[0],
        spans_nm[1],
        spans_nm[2],
    )


def format_particle_size_nm(atoms):
    sx, sy, sz = get_particle_size_nm(atoms)

    return (
        f"{sx:.2f} x "
        f"{sy:.2f} x "
        f"{sz:.2f} nm"
    )


def format_parameter_for_filename(value):
    """
    Format numeric parameters for output filenames without
    losing meaningful decimal precision.

    Examples
    --------
    0.25 -> "0.25"
    0.5  -> "0.5"
    1.0  -> "1.0"
    1.5  -> "1.5"
    """

    text = f"{float(value):.6f}".rstrip("0").rstrip(".")

    if "." not in text:
        text += ".0"

    return text


def get_shell_fraction_percent(
    shell_atoms,
    total_atoms,
):
    if total_atoms == 0:
        return 0.0

    return (
        100.0
        * shell_atoms
        / total_atoms
    )


def run_once():
    print("NanoAmorphizer")
    print(
        "Crystalline nanoparticle generator "
        "with disordered surface shell"
    )
    print()

    # --------------------------------------------------
    # Load CIF
    # --------------------------------------------------

    cif_path = ask_cif_path()

    try:
        atoms = load_cif(cif_path)

    except Exception as exc:
        print()
        print(
            f"Error loading CIF: {exc}"
        )
        return

    print()
    print("CIF loaded successfully.")
    print()

    # --------------------------------------------------
    # Structure summary
    # --------------------------------------------------

    formula = atoms.get_chemical_formula()
    natoms = len(atoms)

    cell = atoms.get_cell()
    lengths = cell.lengths()
    volume = cell.volume

    print("Structure summary")
    print("-----------------")
    print(f"Formula: {formula}")
    print(
        f"Atoms in unit cell: "
        f"{natoms}"
    )
    print(
        f"Cell lengths (A): "
        f"{lengths[0]:.3f}, "
        f"{lengths[1]:.3f}, "
        f"{lengths[2]:.3f}"
    )
    print(
        f"Cell volume (A^3): "
        f"{volume:.3f}"
    )

    # --------------------------------------------------
    # Particle parameters
    # --------------------------------------------------

    print()
    print("Particle parameters")
    print("-------------------")

    try:
        radius = ask_radius()
        validate_radius(radius)

        shell = ask_shell()
        validate_shell(
            shell,
            radius,
        )

        disorder = ask_disorder()
        validate_disorder(
            disorder,
            shell,
        )

        seed = ask_seed()

    except Exception as exc:
        print()
        print(
            f"Input error: {exc}"
        )
        return

    print()
    print("Parameters accepted:")
    print(
        f"Radius (A): "
        f"{radius}"
    )
    print(
        f"Shell thickness (A): "
        f"{shell}"
    )
    print(
        f"Disorder strength (A): "
        f"{disorder}"
    )
    print(
        f"Random seed: "
        f"{seed if seed is not None else 'random'}"
    )

    # --------------------------------------------------
    # Build supercell
    # --------------------------------------------------

    print()
    print("Building supercell...")
    print("---------------------")

    supercell, reps = build_supercell(
        atoms,
        radius,
    )

    nx, ny, nz = reps

    print(
        f"Replication factors: "
        f"{nx} x {ny} x {nz}"
    )
    print(
        f"Total atoms in supercell: "
        f"{len(supercell)}"
    )

    # --------------------------------------------------
    # Crop spherical nanoparticle
    # --------------------------------------------------

    print()
    print("Cropping spherical particle...")
    print("------------------------------")

    particle, particle_center = crop_sphere(
        supercell,
        radius,
    )

    atoms_after_crop = len(particle)

    # --------------------------------------------------
    # Remove isolated surface atoms
    # --------------------------------------------------

    particle, isolated_atoms_removed = (
        remove_isolated_atoms(
            particle,
        )
    )

    # --------------------------------------------------
    # Core-shell analysis
    # --------------------------------------------------

    shell_info = get_shell_info(
        particle,
        radius=radius,
        shell_thickness=shell,
        center=particle_center,
    )

    total_atoms = (
        shell_info["total_atoms"]
    )
    core_atoms = (
        shell_info["core_atoms"]
    )
    shell_atoms = (
        shell_info["shell_atoms"]
    )
    outside_atoms = (
        shell_info["outside_atoms"]
    )

    shell_fraction = (
        get_shell_fraction_percent(
            shell_atoms,
            total_atoms,
        )
    )

    # --------------------------------------------------
    # First-neighbor analysis
    # --------------------------------------------------

    print()
    print("Analyzing local structure...")
    print("----------------------------")

    try:
        neighbor_network = (
            build_neighbor_network(
                particle
            )
        )

        neighbor_summary = (
            get_neighbor_summary(
                particle,
                neighbor_network=neighbor_network,
            )
        )

    except Exception as exc:
        print()
        print(
            "Local structure analysis "
            f"error: {exc}"
        )
        return

    # --------------------------------------------------
    # Particle information
    # --------------------------------------------------

    print()
    print("Nanoparticle summary")
    print("--------------------")
    print(
        f"Atoms after spherical crop: "
        f"{atoms_after_crop}"
    )
    print(
        f"Isolated surface atoms removed: "
        f"{isolated_atoms_removed}"
    )
    print(
        f"Total atoms in particle: "
        f"{total_atoms}"
    )
    print(
        f"Core atoms: "
        f"{core_atoms}"
    )
    print(
        f"Shell atoms: "
        f"{shell_atoms}"
    )
    print(
        f"Outside atoms: "
        f"{outside_atoms}"
    )
    print(
        f"Shell fraction: "
        f"{shell_fraction:.1f} %"
    )
    print(
        f"Estimated nanoparticle size: "
        f"{format_particle_size_nm(particle)}"
    )
    print(
        f"First-neighbor pairs: "
        f"{neighbor_summary['neighbor_pairs']}"
    )
    print(
        "Characteristic nearest-neighbor "
        f"distance: "
        f"{neighbor_summary['characteristic_distance_A']:.3f} A"
    )
    print(
        f"Median coordination: "
        f"{neighbor_summary['median_coordination']:.1f}"
    )

    # --------------------------------------------------
    # Export crystalline nanoparticle
    # --------------------------------------------------

    print()
    print(
        "Exporting crystalline nanoparticle..."
    )
    print(
        "------------------------------------"
    )

    crystalline_filename = (
        f"{formula}_"
        f"R{int(radius)}_"
        f"S{int(shell)}_"
        f"D{format_parameter_for_filename(disorder)}_"
        f"crystalline.xyz"
    )

    crystalline_path = export_xyz(
        particle,
        crystalline_filename,
    )

    print(
        "Crystalline nanoparticle generated."
    )
    print(
        f"Path: {crystalline_path}"
    )

    amorphized_particle = None
    amorphous_path = None
    amorphous_shell_info = None
    amorphization_report = None

    # --------------------------------------------------
    # Amorphization
    # --------------------------------------------------

    print()
    print("Amorphization step")
    print("------------------")

    if shell <= 0:

        print(
            "Shell thickness is 0 A."
        )
        print(
            "No amorphous particle generated."
        )

    elif disorder <= 0:

        print(
            "Disorder strength is 0 A."
        )
        print(
            "No amorphous particle generated."
        )

    else:

        print(
            "Generating amorphous shell..."
        )

        (
            amorphized_particle,
            amorphization_report,
        ) = amorphize_shell(
            particle,
            radius=radius,
            shell_thickness=shell,
            disorder_strength=disorder,
            seed=seed,
            center=particle_center,
            return_report=True,
        )

        amorphous_shell_info = (
            get_shell_info(
                amorphized_particle,
                radius=radius,
                shell_thickness=shell,
                center=particle_center,
            )
        )

        amorphous_filename = (
            f"{formula}_"
            f"R{int(radius)}_"
            f"S{int(shell)}_"
            f"D{format_parameter_for_filename(disorder)}_"
            f"amorphous.xyz"
        )

        amorphous_path = export_xyz(
            amorphized_particle,
            amorphous_filename,
        )

        print(
            "Amorphous nanoparticle generated."
        )
        print(
            f"Path: {amorphous_path}"
        )

    # --------------------------------------------------
    # Generation summary
    # --------------------------------------------------

    print()
    print("Generation summary")
    print("------------------")

    print(
        f"Total atoms: "
        f"{total_atoms}"
    )
    print(
        f"Isolated surface atoms removed: "
        f"{isolated_atoms_removed}"
    )
    print(
        f"Initial core atoms: "
        f"{core_atoms}"
    )
    print(
        f"Initial shell atoms: "
        f"{shell_atoms}"
    )
    print(
        f"Initial outside atoms: "
        f"{outside_atoms}"
    )
    print(
        f"Initial shell fraction: "
        f"{shell_fraction:.1f} %"
    )
    print(
        f"Random seed: "
        f"{seed if seed is not None else 'random'}"
    )
    print(
        f"Crystalline size: "
        f"{format_particle_size_nm(particle)}"
    )
    print(
        f"Nearest-neighbor distance: "
        f"{neighbor_summary['characteristic_distance_A']:.3f} A"
    )
    print(
        f"Crystalline file: "
        f"{crystalline_path}"
    )

    if amorphized_particle is not None:

        print(
            f"Shell atoms selected for disorder: "
            f"{shell_atoms}"
        )
        print(
            f"Atoms actually displaced: "
            f"{amorphization_report['atoms_displaced']}"
        )
        print(
            f"Zero-weight shell atoms: "
            f"{amorphization_report['zero_weight_atoms']}"
        )
        print(
            f"Atoms retained after max attempts: "
            f"{amorphization_report['atoms_retained_after_max_attempts']}"
        )
        print(
            f"Total displacement proposals: "
            f"{amorphization_report['total_proposals']}"
        )
        print(
            f"Rejected proposals: "
            f"{amorphization_report['rejected_proposals']}"
        )
        print(
            f"  overlap: "
            f"{amorphization_report['rejected_overlap']}"
        )
        print(
            f"  parent distance: "
            f"{amorphization_report['rejected_parent_distance']}"
        )
        print(
            f"  forbidden pair type: "
            f"{amorphization_report['rejected_pair_type']}"
        )
        print(
            f"  angle: "
            f"{amorphization_report['rejected_angle']}"
        )
        print(
            f"Proposal acceptance rate: "
            f"{amorphization_report['proposal_acceptance_rate_percent']:.1f} %"
        )
        print(
            f"Shell displacement rate: "
            f"{amorphization_report['shell_displacement_rate_percent']:.1f} %"
        )
        print(
            f"Mean proposals per attempted atom: "
            f"{amorphization_report['mean_proposals_per_attempted_atom']:.2f}"
        )
        print(
            f"Final core-region atoms: "
            f"{amorphous_shell_info['core_atoms']}"
        )
        print(
            f"Final shell-region atoms: "
            f"{amorphous_shell_info['shell_atoms']}"
        )
        print(
            f"Atoms outside nominal radius: "
            f"{amorphous_shell_info['outside_atoms']}"
        )
        print(
            f"Amorphous size: "
            f"{format_particle_size_nm(amorphized_particle)}"
        )
        print(
            f"Amorphous file: "
            f"{amorphous_path}"
        )

    else:

        print(
            "Amorphous file: not generated"
        )


def main():
    while True:

        print()
        run_once()
        print()

        answer = input(
            "Would you like to generate "
            "another nanoparticle? (y/n): "
        ).strip().lower()

        if answer not in [
            "y",
            "yes",
        ]:

            print()
            print(
                "NanoAmorphizer finished."
            )
            break


if __name__ == "__main__":
    main()
