from qiskit import QuantumCircuit
from qiskit.circuit.library import TwoLocal
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager, PassManagerConfig

import numpy as np
import matplotlib.pyplot as plt
from qiskit.circuit.library import RealAmplitudes

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit.visualization import plot_gate_map
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import QiskitRuntimeService

# Define the Hamiltonian

# Create the ansatz
ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz', reps=3)
ansatz = RealAmplitudes(10, reps=1)
"ibm_brisbane"
service = QiskitRuntimeService(channel="ibm_quantum")
backend = service.get_backend("ibm_brisbane")

# pass_manager_config = PassManagerConfig(basis_gates=backend.configuration().basis_gates,
#                                         coupling_map=backend.configuration().coupling_map,
#                                         initial_layout=None)
pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=3)

# Transpile the optimized circuit using the preset pass manager
transpiled_circuit = pass_manager.run(ansatz)
print(ansatz.decompose())
# Get the depth of the transpiled circuit
circuit_depth = transpiled_circuit.depth()
print(f"Optimized Circuit Depth: {circuit_depth}")
