# Mutant builder

Miguel Arbesú -- <miguelarbesu@gmai.com> 

## Description

A script to build mutants from a given PDB file using [Pyrosetta](https://www.pyrosetta.org/). It will replace the existing amino acid at a given position (e.g. 600) with the requested new amino acid in 1-letter code (e.g. "E").

It uses Rossetta's [`FastRelax`](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover) and [`mutate_residue`](https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.toolbox.mutants.html#pyrosetta.toolbox.mutants.mutate_residue) protocols under the hood. 

## Requirements & installation

Your will need Rosetta and PyRosetta licenses to use this script. Both suites are free for academic usage. Your can obtain yours [here](https://els2.comotion.uw.edu/product/pyrosetta).

For installation, the most straitforward is probably via `conda`.

1. [Create a Python environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) . Python 3.6-3.9 is supported. 
2. Follow the instructions [`here`](https://www.pyrosetta.org/downloads#h.c0px19b8kvuw) to add the pyrosetta channel to your profile (license user and pass required)
3. Install `pyrosetta`.

## Usage

```
usage: mutate.py [-h] [--wt_structure WT_STRUCTURE] [--chain CHAIN] [--mutant_position MUTANT_POSITION] [--mutant_aa MUTANT_AA]
                 [--pack_radius PACK_RADIUS]

optional arguments:
  -h, --help            show this help message and exit
  --wt_structure WT_STRUCTURE
                        Path to the wild type structure to mutate
  --chain CHAIN         Chain to mutate. Defaults to A -- the only in AF2 predictions
  --mutant_position MUTANT_POSITION
                        Position to mutate in the canonical sequence
  --mutant_aa MUTANT_AA
                        Amino acid to introduce in 1-letter code
  --pack_radius PACK_RADIUS
                        Radius to determine which neighbor side chains are re-packed
```

## Resources

The script is based in the following exammples and documentation pages.

- https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.toolbox.mutants.html#pyrosetta.toolbox.mutants.mutate_residue
- https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.08-Point-Mutation-Scan.ipynb
- https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/02.02-Working-with-Pose-Residues.ipynb


