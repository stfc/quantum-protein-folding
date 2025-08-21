import os
from itertools import product
import numpy as np
import json
from qiskit.quantum_info import SparsePauliOp


def generate_bitstrings(num_res, num_rot):
    """
    Generates all valid bitstrings with Hamming weight of 1.
    """
    base_strings = ["0" * i + "1" + "0" * (num_rot - i - 1) for i in range(num_rot)]
    bitstrings = ["".join(bits) for bits in product(base_strings, repeat=num_res)]
    return bitstrings


def calculate_bitstring_energies(
    bitstrings: list[str], hamiltonian: SparsePauliOp
) -> dict[str, float]:
    """
    Vectorised computation of energies for a list of bitstrings using a SparsePauliOp Hamiltonian.

    Parameters:
    - bitstrings (List[str]): List of binary strings
    - hamiltonian (SparsePauliOp): Hamiltonian as a qiskit.quantum_info.SparsePauliOp object

    Returns:
    - dict[str, float]: Dictionary mapping each bitstring to its energy
    """
    bit_array = np.array(
        [[int(b) for b in bitstring] for bitstring in bitstrings], dtype=int
    )
    energies = np.zeros(len(bitstrings), dtype=np.complex128)

    for pauli, coeff in zip(hamiltonian.paulis, hamiltonian.coeffs):
        if abs(coeff) < 1e-10:
            continue
        z_mask = np.array([op == "Z" for op in str(pauli)], dtype=bool)
        z_values = 1 - 2 * bit_array[:, z_mask]  # Z on |0> = +1, Z on |1> = -1
        term_expectation = np.prod(z_values, axis=1)
        energies += coeff * term_expectation

    if not np.all(np.isreal(energies)):
        raise ValueError("Energies array contains complex values.")

    energies = np.real(energies)
    return {bitstring: float(energy) for bitstring, energy in zip(bitstrings, energies)}


def get_min_energy_bitstring(filepath: str) -> str:
    """
    Loads a JSON file containing bitstring-energy mappings and returns the bitstring
    with the minimum energy value.

    Parameters:
    - filepath (str): Path to the JSON file

    Returns:
    - str: Bitstring corresponding to the minimum energy
    """

    if not os.path.exists(filepath):
        raise FileNotFoundError(
            "JSON file not found. Generate by executing scripts/calculate_exact_energies.py"
        )

    with open(filepath, "r", encoding='utf-8') as f:
        data = json.load(f)

    return min(data, key=data.get)

def bitstring_to_int(bitstring):
    """Converts a binary bitstring to an integer."""
    return int(bitstring, 2)