from ase.io import read
import os


def load_cif(cif_path):
    if not os.path.isfile(cif_path):
        raise FileNotFoundError(f"File not found: {cif_path}")

    try:
        atoms = read(cif_path)
    except Exception as exc:
        raise ValueError(f"Could not read CIF file: {exc}")

    return atoms