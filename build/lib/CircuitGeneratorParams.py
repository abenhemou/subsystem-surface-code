from typing import Callable, Set, List, Dict, Tuple, Optional
import math
import stim 
from dataclasses import dataclass

def append_anti_basis_error(circuit: stim.Circuit, targets: List[int], p: float, basis: str) -> None:
    if p > 0:
        if basis == "X":
            circuit.append_operation("Z_ERROR", targets, p)
        else:
            circuit.append_operation("X_ERROR", targets, p)

@dataclass
class CircuitGenParameters:
    distance: int
    rounds: int
    after_clifford_depolarization: float = 0
    before_round_data_depolarization: float = 0
    before_measure_flip_probability: float = 0
    after_reset_flip_probability: float = 0
    exclude_other_basis_detectors: bool = False
    show_qubit_coordinates: bool = True

    def append_reset(self, circuit: stim.Circuit, targets: List[int], basis: str = "Z") -> None:
        """Append a reset operation in a desired basis on a given list of qubits."""
        circuit.append_operation("R" + basis, targets)
        append_anti_basis_error(circuit, targets, self.after_reset_flip_probability, basis)

    def append_measure(self, circuit: stim.Circuit, targets: List[int], basis: str = "Z") -> None:
        """Append a measure operation in a desired basis on a given list of qubits."""
        append_anti_basis_error(circuit, targets, self.before_measure_flip_probability, basis)
        circuit.append_operation("M" + basis, targets)

    def append_measure_reset(self, circuit: stim.Circuit, targets: List[int], basis: str = "Z") -> None:
        """Append a measure AND reset operation in a desired basis on a given list of qubits."""
        append_anti_basis_error(circuit, targets, self.before_measure_flip_probability, basis)
        circuit.append_operation("MR" + basis, targets)
        append_anti_basis_error(circuit, targets, self.after_reset_flip_probability, basis)

    def append_begin_round_tick(self, circuit: stim.Circuit, data_qubits: List[int]) -> None:
        """Append time step separation indicator in the circuit, and add depolarising noise on the data wubits if desired before every round."""
        circuit.append_operation("TICK", [])
        if self.before_round_data_depolarization > 0:
            circuit.append_operation("DEPOLARIZE1", data_qubits, self.before_round_data_depolarization)
            
    def append_unitary_1(self, circuit: stim.Circuit, name: str, targets: List[int]) -> None:
        """Append a single_qubit unitary to the circuit applied on target list."""
        circuit.append_operation(name, targets)
        if self.after_clifford_depolarization > 0:
            circuit.append_operation("DEPOLARIZE1", targets, self.after_clifford_depolarization)
        
    def append_unitary_2(self, circuit: stim.Circuit, name: str, targets: List[int]) -> None:
        """Append a two_qubit unitary to the circuit applied on target list."""
        circuit.append_operation(name, targets)
        if self.after_clifford_depolarization > 0:
            circuit.append_operation("DEPOLARIZE2", targets, self.after_clifford_depolarization)

def coord_to_idx(q: complex, distance: int) -> int:
    q = q - math.fmod(q.real, 2)*1j
    return int(q.real + q.imag * (distance + 0.5))