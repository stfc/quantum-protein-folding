import argparse
from path_setup import load_project_root
import csv
import os
from scipy.optimize import dual_annealing

load_project_root()

from quantum_protein_folding.quantum_protein_folding.ising import get_hamiltonian
from quantum_protein_folding.quantum_protein_folding.bitstrings import (
    get_min_energy_bitstring,
    calculate_energies_numpy_sparse,
)
from quantum_protein_folding.quantum_protein_folding.config import (
    EXACT_ENERGY_DATA_DIR,
    SIMULATED_ANNEALING_FILE,
)


class CallTracker:
    def __init__(self):
        self.calls = 0


class FoundGroundState(Exception):
    pass


def parse_args():
    """Command line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Run simulated annealing with specified parameters."
    )

    parser.add_argument(
        "-rot",
        "--num_rot",
        type=int,
        required=True,
        help="Number of rotamers per residue",
    )
    parser.add_argument(
        "-res", "--num_res", type=int, required=True, help="Number of residues"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    tracker = CallTracker()

    hamiltonian = get_hamiltonian(num_rot=args.num_rot, num_res=args.num_res)
    num_bits = args.num_rot * args.num_res
    ground_state_filepath = (
        f"{EXACT_ENERGY_DATA_DIR}/res-{args.num_res}-rot-{args.num_rot}.json"
    )
    ground_state_bitstring = get_min_energy_bitstring(ground_state_filepath)

    def objective_function(x):
        binary_state = [1 if xi > 0.5 else 0 for xi in x]
        binary_state = "".join(str(digit) for digit in binary_state)
        tracker.calls += 1
        energy = calculate_energies_numpy_sparse([binary_state], hamiltonian)[binary_state]
        if binary_state == ground_state_bitstring:
            raise FoundGroundState()
        return energy

    bounds = [(0, 1)] * num_bits
    found = False
    max_restarts = 50
    restarts = 0

    while not found and restarts < max_restarts:
        try:
            _ = dual_annealing(objective_function, bounds)
        except FoundGroundState:
            found = True
        restarts += 1

    total_complexity = tracker.calls * (
        num_bits**2
    )  # each call evaluates a matrix expectation value, with O(N^2) complexity

    result_dict = {
        "num_res": args.num_res,
        "num_rot": args.num_rot,
        "ground_state": ground_state_bitstring,
        "gs_found": found,
        "CPU calls": total_complexity,
    }

    write_file = str(SIMULATED_ANNEALING_FILE)

    file_exists = os.path.isfile(write_file)

    if not file_exists:
        os.makedirs(os.path.dirname(write_file), exist_ok=True)

    with open(write_file, mode="a", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=result_dict.keys())
        if not file_exists:
            writer.writeheader()

        writer.writerow(result_dict)


if __name__ == "__main__":
    main()
