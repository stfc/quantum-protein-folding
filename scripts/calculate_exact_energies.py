import os
import json
import argparse
from path_setup import load_project_root

load_project_root()

from quantum_protein_folding.quantum_protein_folding.processing import (
    generate_bitstrings,
    calculate_bitstring_energies,
)
from quantum_protein_folding.quantum_protein_folding.ising import (
    get_hamiltonian,
    get_q_hamiltonian,
)
from quantum_protein_folding.quantum_protein_folding.config import EXACT_ENERGY_DATA_DIR


def parse_args():

    parser = argparse.ArgumentParser(
        description="Generate bitstrings for rotamer states."
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
    """Calculates the energies for all Hamming weight conserving bitstrings for given
    residue and rotomer values. Results are saved as json."""
    args = parse_args()
    num_res = args.num_res
    num_rot = args.num_rot

    if num_res <= 1 or num_rot <= 1:
        raise ValueError("Number of residues and rotamers must be positive integers.")

    bitstrings = generate_bitstrings(num_res, num_rot)
    hamiltonian = get_hamiltonian(num_rot=num_rot, num_res=num_res)
    q_hamiltonian = get_q_hamiltonian(num_res * num_rot, hamiltonian)
    energies = calculate_bitstring_energies(bitstrings, q_hamiltonian)

    os.makedirs(EXACT_ENERGY_DATA_DIR, exist_ok=True)

    filename = f"res-{num_res}-rot-{num_rot}.json"
    file_path = os.path.join(EXACT_ENERGY_DATA_DIR, filename)

    with open(file_path, "w") as f:
        json.dump(energies, f, indent=4)

    print(f"Bitstring energies saved to {file_path}")


if __name__ == "__main__":
    main()
