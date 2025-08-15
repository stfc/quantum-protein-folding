import argparse
import csv
import os
from qiskit_aer import AerSimulator
from qiskit.primitives import BackendSampler
from qiskit.transpiler import PassManager
from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms.minimum_eigensolvers import QAOA

from path_setup import load_project_root

load_project_root()

from quantum_protein_folding.quantum_protein_folding.ising import (
    get_hamiltonian,
    get_q_hamiltonian,
    get_xy_mixer,
)
from quantum_protein_folding.quantum_protein_folding.qaoa import (
    get_initial_parameters,
    get_symmetry_preserving_initial_state,
)
from quantum_protein_folding.quantum_protein_folding.bitstrings import (
    get_min_energy_bitstring,
    bitstring_to_int,
)
from quantum_protein_folding.quantum_protein_folding.config import (
    EXACT_ENERGY_DATA_DIR,
    MPS_QAOA_FILE,
    SV_QAOA_FILE,
)


class GroundStateTracker:
    """Class to keep track of if/when the ground state is found."""

    def __init__(self, int_ground_state, args):
        self.shots_to_ground_state = 0
        self.iterations = 0
        self.int_ground_state = int_ground_state
        self.num_shots = 0
        self.gs_found = False
        self.args = args

    def callback(self, all_bitstrings):
        if self.int_ground_state in all_bitstrings:
            self.shots_to_ground_state = self.num_shots
            self.gs_found = True

        else:
            self.num_shots += self.args.shots
            self.iterations += 1


def parse_args():
    """Command line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Run quantum protein folding simulation with specified parameters."
    )
    parser.add_argument(
        "-p", "--p", type=int, required=True, help="Number of ansatz layers for QAOA"
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
    parser.add_argument(
        "-alpha",
        "--alpha",
        type=float,
        default=0.02,
        help="Alpha parameter for QAOA or other algorithm",
    )
    parser.add_argument(
        "-s",
        "--shots",
        type=int,
        default=1000,
        help="Number of shots for quantum execution",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=["MPS", "SV"],
        default="MPS",
        help="Simulation method (MPS/SV).",
    )
    return parser.parse_args()


def main():
    ### args
    args = parse_args()
    num_rot = args.num_rot
    num_res = args.num_res
    p = args.p

    ### reference data
    exact_energy_filepath = f"{EXACT_ENERGY_DATA_DIR}/res-{num_res}-rot-{num_rot}.json"
    bitstring_ground_state = get_min_energy_bitstring(filepath=exact_energy_filepath)
    int_ground_state = bitstring_to_int(bitstring_ground_state)

    ### qaoa setup
    hamiltonian = get_hamiltonian(num_rot=num_rot, num_res=num_res)
    q_hamiltonian = get_q_hamiltonian(num_res * num_rot, hamiltonian)
    mixer = get_xy_mixer(num_res=num_res, num_rot=num_rot)
    initial_point = get_initial_parameters(p, cost_bound=0.1, mixer_bound=1.0)
    initial_state = get_symmetry_preserving_initial_state(
        num_res=num_res, num_rot=num_rot
    )

    ### simulator
    match args.method:
        case "MPS":
            backend = "matrix_product_state"
        case "SV":
            backend = "statevector"
    simulator = AerSimulator(method=backend)

    ### sampler
    backend_options = {
        "seed_simulator": 42,
        "shots": args.shots,
        "optimization_level": 3,
        "resilience_level": 3,
    }

    sampler = BackendSampler(
        backend=simulator, options=backend_options, bound_pass_manager=PassManager()
    )

    ### ground state tracking
    tracker = GroundStateTracker(int_ground_state, args)

    ### qaoa
    qaoa = QAOA(
        sampler=sampler,
        optimizer=COBYLA(),
        reps=p,
        initial_state=initial_state,
        mixer=mixer,
        initial_point=initial_point,
        callback=tracker.callback,
        aggregation=args.alpha,
    )
    qaoa.optimizer.set_options(maxiter=10_000)

    result = qaoa.compute_minimum_eigenvalue(q_hamiltonian)

    ### recording
    result_dict = {
        "num_res": args.num_res,
        "num_rot": args.num_rot,
        "alpha": args.alpha,
        "shots": args.shots,
        "p": args.p,
        "ground_state": bitstring_ground_state,
        "gs_found": tracker.gs_found,
        "iterations": tracker.iterations,
        "QPU calls": tracker.shots_to_ground_state if tracker.gs_found else -1,
    }

    write_file = None
    if args.method == "MPS":
        write_file = str(MPS_QAOA_FILE)
    if args.method == "SV":
        write_file = str(SV_QAOA_FILE)

    file_exists = os.path.isfile(write_file)

    with open(write_file, mode="a", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=result_dict.keys())
        if not file_exists:
            writer.writeheader()

        writer.writerow(result_dict)


if __name__ == "__main__":
    main()
