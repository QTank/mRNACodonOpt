from qiskit.circuit.library import RealAmplitudes
from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms.minimum_eigensolvers import SamplingVQE
from qiskit.primitives import Sampler
from qiskit.circuit.library import EfficientSU2
from qiskit.opflow import PauliSumOp


def get_min(qubit_op):
    # set classical optimizer
    optimizer = COBYLA(maxiter=50)

    # set variational ansatz
    ansatz = RealAmplitudes(reps=1)
    # ansatz = EfficientSU2(qubit_op.num_qubits, entanglement='linear')
    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    # initialize VQE using CVaR with alpha = 0.1
    vqe = SamplingVQE(
        Sampler(),
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=0.1,
        callback=store_intermediate_result,
    )
    raw_result = vqe.compute_minimum_eigenvalue(qubit_op)
    # print(vqe.sampler.circuits[0].depth())
    return raw_result.best_measurement['bitstring'], raw_result.best_measurement['value']


def test(qubit_op):
    from qiskit.utils import algorithm_globals
    from qiskit_aer.primitives import Estimator as AerEstimator
    from qiskit_ibm_runtime import Options
    seed = 170
    algorithm_globals.random_seed = seed

    noiseless_estimator = AerEstimator(
        run_options={"seed": seed, "shots": 1024},
        transpile_options={"seed_transpiler": seed},
    )
    # set classical optimizer
    optimizer = COBYLA(maxiter=50)

    # set variational ansatz
    ansatz = RealAmplitudes(reps=1)
    # ansatz = EfficientSU2(qubit_op.num_qubits, entanglement='linear')
    counts = []
    values = []

    backend = "ibmq_qasm_simulator"  # use the simulator
    options = Options()
    options.simulator.seed_simulator = 42
    options.execution.shots = 1000
    options.optimization_level = 0  # no optimization
    options.resilience_level = 0  # no error mitigation

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    # initialize VQE using CVaR with alpha = 0.1
    vqe = SamplingVQE(
        Sampler(backend=backend, Options=options),
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=0.1,
        callback=store_intermediate_result,
    )

    raw_result = vqe.compute_minimum_eigenvalue(PauliSumOp(qubit_op))
    # print(vqe.sampler.circuits[0].depth())
    return raw_result.best_measurement['bitstring'], raw_result.best_measurement['value']
