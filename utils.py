import numpy as np

from typing import Union
from scipy.sparse import csc_matrix, identity
from qulacs.circuit import QuantumCircuitOptimizer as QCO
from qulacs.gate import merge
from qulacs import QuantumCircuit as QC

def circuit_to_unitary(circuit):
    
    # qulacs
    if type(circuit) == QC:
        from qulacs.gate import Identity
        n_qubit = circuit.get_qubit_count()
        qco = QCO()
        circuit_ = circuit.copy()
        for i in range(n_qubit):
            circuit_.add_gate(Identity(i))
        qco.optimize_light(circuit_)
        for ind in range(3, n_qubit+1):
            qco.optimize(circuit_, ind)
        gate_list = [circuit_.get_gate(i) for i in range(circuit_.get_gate_count())]
        merged = merge(gate_list)
        return merged.get_matrix()
    
    else:
        raise NotImplementedError()

def word_to_gate(gateword, backend = "qulacs"):
    phase = np.pi/4 * gateword.count("W")
    
    if backend == "qulacs":
        from qulacs import QuantumCircuit as QC
        circuit = QC(1)
        circuit.add_dense_matrix_gate(0, np.array([[np.exp(1j * phase), 0], [0, np.exp(1j * phase)]]))
        for g in gateword:
            if g=="H":
                circuit.add_H_gate(0)
            elif g=="T":
                #circuit.add_dense_matrix_gate(0, [[1, 0], [0, np.exp(1j * np.pi/4)]])
                circuit.add_T_gate(0)
            elif g=="S":
                #circuit.add_dense_matrix_gate(0, [[1, 0], [0, np.exp(1j * np.pi/2)]])
                circuit.add_S_gate(0)
            elif g=="X":
                circuit.add_X_gate(0)                
        
    elif backend == "qiskit":
        from qiskit import QuantumCircuit
        

        circuit = QuantumCircuit(1, global_phase = phase)
        for g in gateword:
            if g=="H":
                circuit.h(0)
            elif g == "T":
                circuit.t(0)
            elif g=="S":
                circuit.s(0)
            elif g=="X":
                circuit.x(0)
    return circuit

def word_to_unitary(gateword, backend = "qulacs"):
    assert backend in ["qulacs", "qiskit"]
    circuit = word_to_gate(gateword, backend)
    if backend == "qiskit":
        from qiskit.quantum_info import Operator
        unitary = Operator(circuit).data
    elif backend == "qulacs":
        unitary = circuit_to_unitary(circuit, )

    return unitary