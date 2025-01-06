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
from qiskit_transpiler_service.transpiler_service import TranspilerService

# Define the Hamiltonian

# Create the ansatz
ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz', reps=3)
ansatz = RealAmplitudes(15, reps=2)
"ibm_brisbane"
service = QiskitRuntimeService(channel="ibm_quantum")
backend = service.get_backend("ibm_brisbane")

# pass_manager_config = PassManagerConfig(basis_gates=backend.configuration().basis_gates,
#                                         coupling_map=backend.configuration().coupling_map,
#                                         initial_layout=None)
pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=3)

# Transpile the optimized circuit using the preset pass manager
transpiled_circuit = pass_manager.run(ansatz)
circuit = ansatz.decompose()
# Get the depth of the transpiled circuit
circuit_depth = transpiled_circuit.depth()
print(f"Optimized Circuit Depth: {circuit_depth}, the number of gates: {transpiled_circuit.num_nonlocal_gates()}")
print(f"Count: {circuit.size()}, ")

transpiler_ai_false = TranspilerService(
    # Add your code here
    backend_name="ibm_brisbane",
    ai="false",
    optimization_level=3,
)

# service = QiskitRuntimeService(channel="ibm_quantum", token="08eed17dbb30f071687e4fc0dab62658cee49e2180d2eae6ea78e9efe82ef84ff6a932469a26ca5333ffcb76e2337df6f8055c7d8d608c1a6805073cc0f3ae78")
# circuit_ai_false = transpiler_ai_false.run(circuit)
# print(f"Transpiled without AI -> Depth: {circuit_ai_false.depth()}, CNOTs: {circuit_ai_false.num_nonlocal_gates()}")
# circuit_ai_false.draw(fold=-1, output="mpl", scale=0.2)