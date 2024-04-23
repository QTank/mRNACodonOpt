from qiskit.quantum_info import SparsePauliOp
from tunable_parameters import TunableParameters
from codon_optimization import CodonOptimization


def get_op(sequence):
    index = 1
    count = 1
    tune_parameters = TunableParameters(0.3, 0.15 * 6, 1, 130000)
    codonOpt = CodonOptimization(sequence, tune_parameters)
    print(f"{index}/{count} Start to execute VQE on {sequence} ... , requiring {codonOpt.qubit_len} qubits")
    qubit_op = codonOpt.create_qubit_op()

    op_list = []
    for op in qubit_op.primitive:
        op_list.append((op.paulis[0], op.coeffs[0].real))
    return SparsePauliOp.from_list(op_list)

H2_op = SparsePauliOp.from_list(
    [
        ("II", -1.052373245772859),
        ("IZ", 0.39793742484318045),
        ("ZI", -0.39793742484318045),
        ("ZZ", -0.01128010425623538),
        ("XX", 0.18093119978423156),
    ]
)

print(f"Number of qubits: {H2_op.num_qubits}")
H2_op = get_op("FVF")
from qiskit.algorithms.minimum_eigensolvers import NumPyMinimumEigensolver
from qiskit.opflow import PauliSumOp

numpy_solver = NumPyMinimumEigensolver()
result = numpy_solver.compute_minimum_eigenvalue(operator=PauliSumOp(H2_op))
ref_value = result.eigenvalue.real
print(f"Reference value: {ref_value:.5f}")

# define ansatz and optimizer
from qiskit.circuit.library import TwoLocal
from qiskit.algorithms.optimizers import SPSA

iterations = 125
ansatz = TwoLocal(rotation_blocks="ry", entanglement_blocks="cz")
spsa = SPSA(maxiter=iterations)

# define callback
# note: Re-run this cell to restart lists before training
counts = []
values = []


def store_intermediate_result(eval_count, parameters, mean, std):
    counts.append(eval_count)
    values.append(mean)

# define Aer Estimator for noiseless statevector simulation
from qiskit.utils import algorithm_globals
from qiskit_aer.primitives import Estimator as AerEstimator

seed = 170
algorithm_globals.random_seed = seed

noiseless_estimator = AerEstimator(
    run_options={"seed": seed, "shots": 1024},
    transpile_options={"seed_transpiler": seed},
)

# instantiate and run VQE
from qiskit.algorithms.minimum_eigensolvers import VQE

vqe = VQE(
    noiseless_estimator, ansatz, optimizer=spsa, callback=store_intermediate_result
)
result = vqe.compute_minimum_eigenvalue(operator=H2_op)

print(f"VQE on Aer qasm simulator (no noise): {result.eigenvalue.real:.5f}")
print(
    f"Delta from reference energy value is {(result.eigenvalue.real - ref_value):.5f}"
)


