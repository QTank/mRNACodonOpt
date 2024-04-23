from qiskit.circuit.library import EfficientSU2

ansatz = EfficientSU2(17, entanglement='linear', reps=3)

print(ansatz.decompose().num_nonlocal_gates())
print(ansatz.decompose().depth())