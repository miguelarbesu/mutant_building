#!/home/marbesu/.conda/envs/pyrosetta/bin/python
# -*- coding: utf-8 -*-

from pathlib import Path
from pyrosetta.toolbox.mutants import mutate_residue
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--wt_structure", type=str, help="Path to the wild type structure to mutate",
    default="data/AF2-trimmed/AF-P15056-F1-model_v3_trimmed.pdb"
)
parser.add_argument(
    "--chain",
    type=str,
    default="A",
    required=False,
    help="Chain to mutate. Defaults to A -- the only in AF2 predictions",
)
parser.add_argument(
    "--mutant_position", type=int, help="Position to mutate in the canonical sequence",
    default=600
)
parser.add_argument(
    "--mutant_aa", type=str, help="Amino acid to introduce in 1-letter code", 
    default="E"
)
parser.add_argument(
    "--pack_radius",
    type=float,
    default=8.0,
    required=False,
    help="Radius to determine which neighbor side chains are re-packed",
)

args = parser.parse_args()

# initialize rosetta
init()

# load wilt type pdb
wt_structure = Path(args.wt_structure)
basePose = pose_from_pdb(str(wt_structure))

relaxed_wt_structure = wt_structure.parent / (wt_structure.stem + "_relaxed.pdb")

# load default full-atom scoring function
scorefxn = get_fa_scorefxn()
# set up relaxation routine
# FastRelax: keep backbone and repack side chains.
# https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover
relax = rosetta.protocols.relax.FastRelax()
relax.set_scorefxn(scorefxn)
relax.constrain_relax_to_start_coords(True)
# limit number of iterations for less resource usage
# relax.max_iter(100)
relax.apply(basePose)
# write the relaxed structure
basePose.dump_pdb(str(relaxed_wt_structure))

# Load the previously-relaxed pose and make a copy to mutate
relaxPose = pose_from_pdb(str(relaxed_wt_structure))


def define_mutation(pose, pdb_mutant_position, mutant_aa, chain="A"):
    # translate canonical numbering (PDB) to pose (shorter)
    pose_mutant_position = pose.pdb_info().pdb2pose(chain, 600)
    wildtype_aa = pose.residue(524).name1()
    if wildtype_aa == mutant_aa:
        print("Mutated amino acid is the same as wild type, please check")
    mutation = wildtype_aa + str(pdb_mutant_position) + mutant_aa
    print(f"Introducing mutation {mutation}")
    return pose_mutant_position, mutation


# define mutation position in the pose  (different from PDB!) and out file
pose_mutant_position, mutation = define_mutation(
    relaxPose, args.mutant_position, args.mutant_aa
)
mutated_structure = relaxed_wt_structure.parent / (
    relaxed_wt_structure.stem + f"_{mutation}.pdb"
)

# make a copy of the relaxed pose and introduce the mutation
mutatedPose = relaxPose.clone()
mutate_residue(
    mutatedPose,
    mutant_position=pose_mutant_position,
    mutant_aa=args.mutant_aa,
    pack_radius=args.pack_radius,
    pack_scorefxn=scorefxn,
)
# write out mutated structure
mutatedPose.dump_pdb(str(mutated_structure))

# report
print(
    f"Original structure energy: {scorefxn(basePose):.3f}",
    f"Relaxed structure energy: {scorefxn(relaxPose):.3f}",
    f"Mutated structure energy: {scorefxn(mutatedPose):.3f}",
    f"RMSD wild type-mutant: {all_atom_rmsd(relaxPose, mutatedPose):.3f}",
    sep="\n",
)
