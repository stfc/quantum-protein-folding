import os
import json
import argparse
from path_setup import load_project_root

load_project_root()

from quantum_protein_folding.quantum_protein_folding.bitstrings import (
    calculate_ground_state_branch_and_bound,
)
from quantum_protein_folding.quantum_protein_folding.ising import get_hamiltonian
from quantum_protein_folding.quantum_protein_folding.config import EXACT_ENERGY_DATA_DIR


def parse_args():

    parser = argparse.ArgumentParser(
        description="Find the ground state energy and bitstring via branch-and-bound."
    )
    parser.add_argument(
        "-res", "--num_res", type=int, required=True, help="Number of residues."
    )
    parser.add_argument(
        "-rot",
        "--num_rot",
        type=int,
        required=True,
        help="Number of rotamers per residue.",
    )
    return parser.parse_args()


def main():
    """Finds the ground state bitstring and energy for given residue and rotamer values
    using branch-and-bound. Result is saved as json."""
    args = parse_args()
    num_res = args.num_res
    num_rot = args.num_rot

    if num_res <= 1 or num_rot <= 1:
        raise ValueError("Number of residues and rotamers must be greater than 1.")

    hamiltonian = get_hamiltonian(num_rot=num_rot, num_res=num_res)
    gs_bitstring, gs_energy = calculate_ground_state_branch_and_bound(
        num_res, num_rot, hamiltonian
    )

    os.makedirs(EXACT_ENERGY_DATA_DIR, exist_ok=True)

    filename = f"res-{num_res}-rot-{num_rot}.json"
    file_path = os.path.join(EXACT_ENERGY_DATA_DIR, filename)

    with open(file_path, "w", encoding="utf-8") as f:
        json.dump({gs_bitstring: gs_energy}, f, indent=4)

    print(f"Ground state saved to {file_path}")
    print(f"  Bitstring: {gs_bitstring}")
    print(f"  Energy:    {gs_energy:.6f}")


if __name__ == "__main__":
    main()
