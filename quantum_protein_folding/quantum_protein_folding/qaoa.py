from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator


### initial params ###
def get_initial_parameters(p, cost_bound=0.1, mixer_bound=1.0):
    """Create random initial parameter vector respecting the cost and mixer unitary bounds."""
    init_point_cost = np.random.uniform(-cost_bound, cost_bound, p)
    init_point_mixer = np.random.uniform(-mixer_bound, mixer_bound, p)
    initial_point = np.zeros(2 * p)
    initial_point[0::2] = init_point_cost
    initial_point[1::2] = init_point_mixer
    return initial_point


### initial state prep ###
def r_gate(theta=np.pi / 4):
    """
    Constructs the R-gate as defined in:

    B. T. Gard et al., "Efficient symmetry-preserving state preparation circuits
    for the variational quantum eigensolver algorithm," npj Quantum Information,
    vol. 6, article 10, 2020.
    """

    qc = QuantumCircuit(1)
    qc.ry(theta * np.pi / 2, 0)
    qc.rz(np.pi, 0)
    op = Operator(qc)
    return op


def a_gate(theta=np.pi / 4):
    """
    Constructs the A-gate as defined in:

    B. T. Gard et al., "Efficient symmetry-preserving state preparation circuits
    for the variational quantum eigensolver algorithm," npj Quantum Information,
    vol. 6, article 10, 2020.
    """
    qc = QuantumCircuit(2)
    qc.cx(1, 0)
    rgate = r_gate(theta=theta)
    rgate_adj = r_gate().adjoint()
    # apply rgate to qubit 0
    qc.unitary(rgate_adj, [1], label="R")
    qc.cx(0, 1)
    qc.unitary(rgate, [1], label="R")
    qc.cx(1, 0)
    return qc


def create_product_state(base_circuit, n):
    """Constructs a product circuit of identical sub circuits."""
    num_qubits = base_circuit.num_qubits
    product_circuit = QuantumCircuit(num_qubits * n)
    for i in range(n):
        base_circuit_copy = deepcopy(base_circuit)
        product_circuit.compose(
            base_circuit_copy,
            qubits=range(i * num_qubits, (i + 1) * num_qubits),
            inplace=True,
        )

    return product_circuit


def get_symmetry_preserving_initial_state(num_res, num_rot, theta=np.pi / 4):
    """Constructs an initial state respecting the required Hamming weight condition
    of the rotomer qubit registeries."""
    if num_rot < 2:
        raise ValueError("num_rot must be at least 2.")

    qc = QuantumCircuit(num_rot)
    qc.x(num_rot // 2)
    agate = a_gate(theta=theta)

    for i in range(num_rot - 1):
        qc.compose(agate, [i, i + 1], inplace=True)

    init_state = create_product_state(qc, num_res)
    return init_state
