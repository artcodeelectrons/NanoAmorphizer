import numpy as np


def build_supercell(atoms, radius):

    cell = atoms.get_cell()
    lengths = cell.lengths()

    target_size = 2 * radius

    nx = int(np.ceil(target_size / lengths[0]))
    ny = int(np.ceil(target_size / lengths[1]))
    nz = int(np.ceil(target_size / lengths[2]))

    supercell = atoms.repeat((nx, ny, nz))

    return supercell, (nx, ny, nz)