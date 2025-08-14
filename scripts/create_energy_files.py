"""This script creates interaction energy files for protein folding side chain optimisation using PyRosetta."""

from pathlib import Path
from pyrosetta import init, pose_from_pdb, get_score_function
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.rotamer_set import RotamerSets
from pyrosetta.rosetta.core.pack.interaction_graph import InteractionGraphFactory
from pyrosetta.rosetta.protocols.relax import FastRelax
from path_setup import load_project_root

load_project_root()

import pandas as pd
import argparse
import pyrosetta
import contextlib
import io


from quantum_protein_folding.quantum_protein_folding.config import (
    PDB_DIR,
    PYROSETTA_ENERGY_FILES_DIR,
)


def parse_args():
    """Command line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Create energy files for protein folding."
    )
    parser.add_argument(
        "-p", "--pdb", required=True, help="PDB filename (without .pdb extension)"
    )
    parser.add_argument(
        "-r",
        "--num_rot",
        type=int,
        default=10,
        help="Number of rotamers to record per residue",
    )
    parser.add_argument(
        "-i",
        "--index",
        type=int,
        default=10,
        help="Starting index for rotamer selection",
    )
    return parser.parse_args()


def prepare_pose(pdb_path):
    """Prepare the pose from a PDB file. The path assumed to be relative to the PDB_DIR. pdb extension is not required."""
    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        pose = pose_from_pdb(str(pdb_path))
        sfxn = get_score_function(True)
        relax = FastRelax()
        relax.set_scorefxn(sfxn)
        relax.apply(pose)
    return pose, sfxn


def build_rotamers(pose, sfxn):
    task = TaskFactory.create_packer_task(pose)
    pose.update_residue_neighbors()
    sfxn.setup_for_packing(pose, task.designing_residues(), task.designing_residues())
    graph = pyrosetta.rosetta.core.pack.create_packer_graph(pose, sfxn, task)

    rotsets = RotamerSets()
    rotsets.set_task(task)
    rotsets.build_rotamers(pose, sfxn, graph)
    rotsets.prepare_sets_for_packing(pose, sfxn)

    ig = InteractionGraphFactory.create_interaction_graph(
        task, rotsets, pose, sfxn, graph
    )
    rotsets.compute_energies(pose, sfxn, graph, ig, 1)

    return rotsets, ig


def extract_two_body_data(pose, rotsets, ig, num_rot, index):
    """Return two-body rotomer interaction energies."""
    residue_count = pose.total_residue()
    data_list = []

    for res_i in range(1, residue_count):
        rot_set_i = rotsets.rotamer_set_for_residue(res_i)
        rot_set_j = rotsets.rotamer_set_for_residue(res_i + 1)
        if rot_set_i is None or rot_set_j is None:
            continue

        molten_i = rotsets.resid_2_moltenres(res_i)
        molten_j = rotsets.resid_2_moltenres(res_i + 1)

        if not ig.find_edge(molten_i, molten_j):
            continue

        for rot_i in range(index, index + num_rot):
            for rot_j in range(index, index + num_rot):
                energy = ig.get_two_body_energy_for_edge(
                    molten_i, molten_j, rot_i, rot_j
                )
                data_list.append(
                    {
                        "res i": res_i,
                        "res j": res_i + 1,
                        "rot A_i": rot_i,
                        "rot B_j": rot_j,
                        "E_ij": energy,
                    }
                )

    df = pd.DataFrame(data_list)
    return df, residue_count


def extract_one_body_data(rotsets, ig, num_rot, index, residue_count):
    """Return on-site rotomer interaction energies."""
    data_list = []

    for res_i in range(1, residue_count + 1):
        rot_set_i = rotsets.rotamer_set_for_residue(res_i)
        if rot_set_i is None:
            continue

        molten_i = rotsets.resid_2_moltenres(res_i)

        for rot_i in range(index, index + num_rot):
            energy = ig.get_one_body_energy_for_node_state(molten_i, rot_i)
            data_list.append({"res i": res_i, "rot A_i": rot_i, "E_ii": energy})

    return pd.DataFrame(data_list)


def main():
    args = parse_args()
    if not args.pdb.endswith(".pdb"):
        args.pdb += ".pdb"
    pdb_path = Path(PDB_DIR) / f"{args.pdb}"

    if not pdb_path.exists():
        raise FileNotFoundError(f"{pdb_path} does not exist.")

    # Suppress PyRosetta init output
    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        init()

    pose, sfxn = prepare_pose(pdb_path)
    rotsets, ig = build_rotamers(pose, sfxn)

    df_two_body, residue_count = extract_two_body_data(
        pose, rotsets, ig, args.num_rot, args.index
    )
    df_one_body = extract_one_body_data(
        rotsets, ig, args.num_rot, args.index, residue_count
    )

    num_res = residue_count
    df_two_body.to_csv(
        f"{PYROSETTA_ENERGY_FILES_DIR}/{args.num_rot}rot_{num_res}res_two_body_terms.csv",
        index=False,
    )
    df_one_body.to_csv(
        f"{PYROSETTA_ENERGY_FILES_DIR}/{args.num_rot}rot_{num_res}res_one_body_terms.csv",
        index=False,
    )

    df_two_body.assign(abs_E=df_two_body["E_ij"].abs()).nlargest(2, "abs_E").drop(
        columns=["abs_E"]
    ).to_csv(f"{PYROSETTA_ENERGY_FILES_DIR}/top_two_body_terms.csv", index=False)

    print("Energy data saved successfully.")


if __name__ == "__main__":
    main()
