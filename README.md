# NanoAmorphizer

**Rapid geometric generation of finite nanoparticle models for scientific visualization.**

NanoAmorphizer is a Python tool for generating crystalline nanoparticles with controllable structurally disordered surface shells

It is designed as a fast and accessible scientific visualization tool for constructing finite nanoparticle models from crystallographic information.

Generated structures are exported as `.xyz` files and can be visualized using software such as **VESTA**, **OVITO**, or other atomistic viewers.

## Main features

* Load crystal structures from CIF files.
* Generate finite spherical crystalline nanoparticles.
* Define nanoparticle radius.
* Generate an optional structurally disordered surface shell.
* Control shell thickness and disorder strength.
* Preserve local geometric coherence through distance, neighbor, and angular constraints.
* Remove isolated atoms with coordination number CN = 0.
* Use reproducible stochastic disorder through random-seed control.
* Export crystalline and modified structures as XYZ files.
* Display structural and generation information through a command-line interface.

## Installation

NanoAmorphizer requires Python 3.8 or newer.

Clone or download the repository and install the required dependencies:

```bash
pip install -r requirements.txt
```

The principal dependencies are:

```text
numpy
ase
```

## Basic usage

From the NanoAmorphizer project directory:

```bash
cd src
python cli.py
```

The program will request the CIF structure and the parameters required for nanoparticle generation.

The standard workflow consists of:

1. Loading a CIF crystal structure.
2. Reviewing the crystallographic information detected by NanoAmorphizer.
3. Selecting the nanoparticle radius.
4. Selecting the shell thickness.
5. Selecting the disorder strength.
6. Generating the crystalline nanoparticle and its structurally disordered counterpart.
7. Exporting the resulting structures as XYZ files.

## Input

NanoAmorphizer uses crystallographic structures provided in `.cif` format.

The CIF file supplies the parent crystal structure used to construct the finite nanoparticle model.

## Output

Generated atomic structures are exported in `.xyz` format.

Depending on the selected workflow, NanoAmorphizer can generate:

* A crystalline spherical nanoparticle;
* A nanoparticle containing a crystalline core and a structurally disordered surface shell.

The resulting XYZ files can be opened directly in compatible atomistic visualization software.

## Model parameters

**Nanoparticle radius**
Controls the approximate size of the generated spherical nanoparticle.

**Shell thickness**
Defines the thickness of the surface region selected for structural disorder.

**Disorder strength**
Controls the magnitude of the positional perturbation applied within the selected shell.

**Random seed**
Allows stochastic shell generation to be reproduced.

## Geometric constraints

NanoAmorphizer v1.0.0 applies geometric constraints during shell generation to preserve local structural coherence while introducing controlled disorder.

The current implementation includes:

* Minimum interatomic-distance control;
* Preservation of parent first-neighbor relationships;
* Restrictions on neighbor chemical types;
* Local angular-deviation control;
* Rejection of displacement proposals that violate the selected geometric criteria.

The constants used by this release are documented in:

```text
docs/CONSTRAINT_CONSTANTS.md
```

## Visualization

NanoAmorphizer does not provide its own atomistic visualization environment.

Generated XYZ structures can be visualized using external software such as:

* VESTA
* OVITO
* Other compatible atomic-structure viewers

## Model scope

NanoAmorphizer is intended for rapid geometric structure generation and scientific visualization.

The generated structures are qualitative geometric models. NanoAmorphizer does not perform first-principles calculations, energy minimization, molecular dynamics, or thermodynamic stability calculations.

## Documentation

Additional technical documentation is available in the `docs/` directory.

## Version

Current public release:

**NanoAmorphizer v1.0.0**

This is the first version prepared for public distribution.

## Citation

If you use NanoAmorphizer in academic or scientific work, please cite the software using the metadata provided in:

```text
CITATION.cff
```

## License

NanoAmorphizer is released under the **GNU General Public License v3.0 (GPL-3.0)**.

Copyright (C) 2026
Miriam Carolina Mendoza-Ramirez

You are free to use, modify, and redistribute this software under the terms of the GPL-3.0 license.

See the `LICENSE` file for the complete license text.

## Author

**Miriam Carolina Mendoza-Ramirez**

NanoAmorphizer
2026
