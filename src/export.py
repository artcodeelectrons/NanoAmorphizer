import os
from ase.io import write


def export_xyz(atoms, output_filename, output_dir="../outputs"):
    """
    Export ASE Atoms object to XYZ format inside the outputs folder.
    If the filename already exists, create a new one with an incremental suffix.

    Parameters
    ----------
    atoms : ASE Atoms
        Structure to export.
    output_filename : str
        Desired xyz filename, for example:
        "O3SrTi_R50_S5_D1.0_crystalline.xyz"
    output_dir : str
        Output folder path relative to src.

    Returns
    -------
    output_path : str
        Final path used to save the file.
    """

    os.makedirs(output_dir, exist_ok=True)

    name, ext = os.path.splitext(output_filename)

    final_filename = output_filename
    output_path = os.path.join(output_dir, final_filename)

    counter = 1
    while os.path.exists(output_path):
        final_filename = f"{name}_{counter:02d}{ext}"
        output_path = os.path.join(output_dir, final_filename)
        counter += 1

    write(output_path, atoms)

    return output_path