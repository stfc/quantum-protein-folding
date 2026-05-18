import os
from itertools import product
import numpy as np
import json
from qiskit.quantum_info import SparsePauliOp


def generate_bitstrings(num_res, num_rot):
    """Generates all valid bitstrings with Hamming weight of 1 per residue block."""
    base_strings = ["0" * i + "1" + "0" * (num_rot - i - 1) for i in range(num_rot)]
    bitstrings = ["".join(bits) for bits in product(base_strings, repeat=num_res)]
    return bitstrings


def calculate_energies_pauli(
    bitstrings: list[str], hamiltonian: SparsePauliOp
) -> dict[str, float]:
    """Computes energies for all bitstrings via SparsePauliOp expectation values (Qiskit)."""
    bit_array = np.array(
        [[int(b) for b in bitstring] for bitstring in bitstrings], dtype=int
    )
    energies = np.zeros(len(bitstrings), dtype=np.complex128)

    for pauli, coeff in zip(hamiltonian.paulis, hamiltonian.coeffs):
        if abs(coeff) < 1e-10:
            continue
        z_mask = np.array([op == "Z" for op in pauli.to_label()], dtype=bool)
        z_values = 1 - 2 * bit_array[:, z_mask]
        term_expectation = np.prod(z_values, axis=1)
        energies += coeff * term_expectation

    if not np.all(np.isreal(energies)):
        raise ValueError("Energies array contains complex values.")

    energies = np.real(energies)
    return {bitstring: float(energy) for bitstring, energy in zip(bitstrings, energies)}


def calculate_energies_numpy_sparse(
    bitstrings: list[str], H: np.ndarray
) -> dict[str, float]:
    """Computes energies for all bitstrings directly from the classical Hamiltonian matrix.

    For each bitstring x with spins s = 1 − 2x (s_i ∈ {-1, +1}):
        E(x) = s·diag(H) + s·triu(H,1)·s

    Equivalent to calculate_energies_pauli but avoids Qiskit entirely via a
    single vectorised numpy matmul. Prefer this over calculate_energies_pauli
    unless a SparsePauliOp is already available.
    """
    bit_array = np.array([[int(b) for b in bs] for bs in bitstrings], dtype=np.float64)
    s = 1.0 - 2.0 * bit_array
    H_diag = np.diag(H)
    H_upper = np.triu(H, 1)
    energies = (s * H_diag).sum(axis=1) + (s @ H_upper * s).sum(axis=1)
    return {bs: float(e) for bs, e in zip(bitstrings, energies)}


def calculate_ground_state_branch_and_bound(
    num_res: int, num_rot: int, H: np.ndarray
) -> tuple[str, float]:
    """Finds the ground state via branch-and-bound over the one-hot configuration space.

    Residues are assigned one at a time. At each node a lower bound is computed
    from (a) the minimum achievable one-body + cross-term for each remaining
    residue given all choices made so far, plus (b) a precomputed suffix sum of
    per-pair minima for all not-yet-assigned residue pairs. Branches whose bound
    cannot beat the current best are pruned.

    Prefer this over calculate_energies_numpy_sparse when only the minimum is needed,
    especially for large problems where full enumeration is infeasible.
    """
    H_diag = np.diag(H)
    H_upper = np.triu(H, 1)

    # h1[i, r]: one-body energy from block i when residue i picks rotamer r (s_r=-1, others=+1)
    h1 = np.zeros((num_res, num_rot))
    for i in range(num_res):
        block = H_diag[i * num_rot : (i + 1) * num_rot]
        s = np.ones(num_rot)
        for r in range(num_rot):
            s[:] = 1.0
            s[r] = -1.0
            h1[i, r] = (block * s).sum()

    # h2[i, r_i, j, r_j]: two-body contribution from cross terms between blocks i and j (i < j)
    h2 = np.zeros((num_res, num_rot, num_res, num_rot))
    for i in range(num_res):
        for j in range(i + 1, num_res):
            block = H_upper[i * num_rot : (i + 1) * num_rot,
                            j * num_rot : (j + 1) * num_rot]
            s_i = np.ones(num_rot)
            s_j = np.ones(num_rot)
            for r_i in range(num_rot):
                s_i[:] = 1.0
                s_i[r_i] = -1.0
                for r_j in range(num_rot):
                    s_j[:] = 1.0
                    s_j[r_j] = -1.0
                    h2[i, r_i, j, r_j] = (block * np.outer(s_i, s_j)).sum()

    # pair_lb[k]: lower bound on pair interactions among residues k, k+1, ..., num_res-1
    pair_min = np.zeros((num_res, num_res))
    for i in range(num_res):
        for j in range(i + 1, num_res):
            pair_min[i, j] = h2[i, :, j, :].min()

    pair_lb = np.zeros(num_res + 1)
    for k in range(num_res - 1, -1, -1):
        pair_lb[k] = pair_lb[k + 1] + pair_min[k, k + 1 :].sum()

    best_energy = [float("inf")]
    best_choices: list[list[int]] = [[]]

    def branch(depth: int, choices: list[int], partial_e: float) -> None:
        if depth == num_res:
            if partial_e < best_energy[0]:
                best_energy[0] = partial_e
                best_choices[0] = list(choices)
            return

        for r in range(num_rot):
            delta = h1[depth, r]
            for prev in range(depth):
                delta += h2[prev, choices[prev], depth, r]
            new_e = partial_e + delta

            remaining_lb = 0.0
            for d in range(depth + 1, num_res):
                costs = h1[d, :].copy()
                for prev in range(depth):
                    costs += h2[prev, choices[prev], d, :]
                costs += h2[depth, r, d, :]
                remaining_lb += costs.min()
            remaining_lb += pair_lb[depth + 1]

            if new_e + remaining_lb < best_energy[0]:
                choices.append(r)
                branch(depth + 1, choices, new_e)
                choices.pop()

    branch(0, [], 0.0)

    bs = "".join(
        "0" * r + "1" + "0" * (num_rot - r - 1) for r in best_choices[0]
    )
    return bs, best_energy[0]


def get_min_energy_bitstring(filepath: str) -> str:
    """Returns the minimum-energy bitstring from a JSON file of {bitstring: energy} mappings."""
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
