import os
import sys
from copy import deepcopy
import numpy as np
import pandas as pd
from qiskit.quantum_info import Pauli, SparsePauliOp
from .config import (
    PYROSETTA_ENERGY_IMMUTABLE_FILES_DIR,
    PYROSETTA_ENERGY_FILES_DIR,
)


def get_hamiltonian(num_rot: int, num_res: int) -> np.ndarray:
    """
    Generates or loads the Hamiltonian matrix for a given number of rotamers and residues.

    Parameters:
    - num_rot (int): Number of rotamers
    - num_res (int): Number of residues

    Returns:
    - np.ndarray: Hamiltonian matrix
    """
    one_body_file = os.path.join(
        PYROSETTA_ENERGY_IMMUTABLE_FILES_DIR,
        f"{num_rot}rot_{num_res}res_one_body_terms.csv",
    )
    two_body_file = os.path.join(
        PYROSETTA_ENERGY_IMMUTABLE_FILES_DIR,
        f"{num_rot}rot_{num_res}res_two_body_terms.csv",
    )

    # Check for missing files
    missing_files = []
    if not os.path.exists(one_body_file):
        missing_files.append("one_body_file")
    if not os.path.exists(two_body_file):
        missing_files.append("two_body_file")

    if missing_files:
        print(
            f"Missing energy files ({', '.join(missing_files)}) "
            f"for num_rot={num_rot}, num_res={num_res}. "
            "Please generate them by running create_energy_files.py"
        )
        sys.exit()

    # Load one-body terms
    one_body_df = pd.read_csv(one_body_file)
    q = one_body_df["E_ii"].values
    num_qubits = len(q)

    # Load two-body terms
    two_body_df = pd.read_csv(two_body_file)
    values = two_body_df["E_ij"].values
    Q = np.zeros((num_qubits, num_qubits))

    idx = 0
    for j in range(0, num_qubits - num_rot, num_rot):
        for i in range(j, j + num_rot):
            for offset in range(num_rot):
                k = j + num_rot + offset
                Q[i, k] = Q[k, i] = deepcopy(values[idx])
                idx += 1

    # Construct Hamiltonian matrix
    H = np.zeros((num_qubits, num_qubits))
    for i in range(num_qubits):
        H[i, i] = -0.5 * q[i] - 0.25 * np.sum(Q[i, np.arange(num_qubits) != i])
        for j in range(num_qubits):
            if i != j:
                H[i, j] = 0.25 * Q[i, j]

    return H


def get_q_hamiltonian(num_qubits: int, h_matrix: np.ndarray) -> SparsePauliOp:
    """
    Converts a classical Hamiltonian matrix into a quantum Hamiltonian
    represented as a SparsePauliOp using Ising encoding.

    Parameters:
    - num_qubits (int): Number of qubits
    - h_matrix (np.ndarray): Classical Hamiltonian matrix

    Returns:
    - SparsePauliOp: Hamiltonian as a qiskit.quantum_info.SparsePauliOp object
    """
    pauli_list = []
    coeffs = []

    for i in range(num_qubits):
        for j in range(i, num_qubits):
            if h_matrix[i][j] != 0:
                z_term = ["I"] * num_qubits
                z_term[i] = "Z"
                if i != j:
                    z_term[j] = "Z"
                pauli_list.append("".join(z_term))
                coeffs.append(h_matrix[i][j])

    # Add identity term with zero coefficient
    pauli_list.append("I" * num_qubits)
    coeffs.append(0.0)

    return SparsePauliOp(pauli_list, coeffs)


def get_xy_mixer(num_res, num_rot, transverse_field=1):
    if num_rot < 2:
        raise ValueError("num_rot must be at least 2.")
    num_qubits = num_rot * num_res
    hamiltonian = SparsePauliOp(Pauli("I" * num_qubits), coeffs=[0])

    for i in range(0, num_qubits - num_rot + 1, num_rot):
        for j in range(num_rot):
            for k in range(j + 1, num_rot):
                xx_term = ["I"] * num_qubits
                yy_term = ["I"] * num_qubits

                xx_term[i + j] = "X"
                xx_term[i + k] = "X"
                yy_term[i + j] = "Y"
                yy_term[i + k] = "Y"

                xx_op = SparsePauliOp(Pauli("".join(xx_term)), coeffs=[1 / 2])
                yy_op = SparsePauliOp(Pauli("".join(yy_term)), coeffs=[1 / 2])

                hamiltonian += xx_op + yy_op

    hamiltonian *= transverse_field
    return hamiltonian
