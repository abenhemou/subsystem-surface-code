import numpy as np
import stim 
import pymatching 
import sinter 
import math

from scipy.optimize import curve_fit
from CircuitGeneratorParams import *
from generate_circuit import *

def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    """Method to sample from a given stim circuit, build its detector graph
       and run minimum-weight perfect matching for a given number of shots."""
    
    num_detectors = circuit.num_detectors
    num_observables = circuit.num_observables

    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Extract decoder configuration data from the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=True)

    # Run the decoder : this is what calls pymatching 
    predictions = sinter.predict_observables(
        dem=detector_error_model,
        dets=detection_events,
        decoder='pymatching',
    ) 

    # Count the mistakes.
    num_errors = 0
    for actual_flip, predicted_flip in zip(observable_flips, predictions):
        if not np.array_equal(actual_flip, predicted_flip):
            num_errors += 1
            
    return num_errors

def ler_sample_from_subsystem_surface_code(distance: int, num_shots: int, p_pauli: float = 0, rounds:int = 2):
    """Sample from subsystem surface code circuit."""

    circuit = generate_rotated_subsystem_surface_code(CircuitGenParameters(distance = distance,
                                                                           rounds = rounds,
                                                                           after_clifford_depolarization = p_pauli,
                                                                           before_round_data_depolarization = p_pauli,
                                                                           before_measure_flip_probability = p_pauli,
                                                                           after_reset_flip_probability = p_pauli,
                                                                           exclude_other_basis_detectors = False,
                                                                           show_qubit_coordinates=False),
                                                      SubsystemSurfaceCode(distance=distance),
                                                      experiment_basis='X')
        
    num_errors = count_logical_errors(circuit, num_shots)
    return num_errors / num_shots