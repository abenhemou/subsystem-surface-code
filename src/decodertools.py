import random
import numpy as np
import stim 
import sinter
import pymatching 
from typing import Callable, Set, List, Dict, Tuple, Optional
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 

from .CircuitGeneratorParams import *
from .generate_circuit import *

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

def ler_sample_from_subsystem_surface_code(distance: int,
                                            num_shots: int,
                                            noise_params: List[float],
                                            rounds:int = 2,
                                            ):
    """Sample from subsystem surface code circuit."""

    circuit = generate_rotated_subsystem_surface_code(CircuitGenParameters(distance = distance,
                                                                           rounds = rounds,
                                                                           after_clifford_depolarization = noise_params[0],
                                                                           before_round_data_depolarization = noise_params[1],
                                                                           before_measure_flip_probability = noise_params[2],
                                                                           after_reset_flip_probability = noise_params[3],
                                                                           exclude_other_basis_detectors = False,
                                                                           show_qubit_coordinates=False),
                                                      SubsystemSurfaceCode(distance=distance),
                                                      experiment_basis='X')
        
    num_errors = count_logical_errors(circuit, num_shots)
    return num_errors / num_shots


#-------------------
# Fitting functions

def get_fit_params(xdata, ydata, params_0=None, ftol: float = 1e-5, maxfev: int = 2000):
    """Get fitting params."""

    # Curve fit.
    bounds = [min(xdata[1]), max(xdata[1])]

    if params_0 is not None and params_0[0] not in bounds:
        params_0[0] = (bounds[0] + bounds[1]) / 2

    params_opt, _ = curve_fit(fit_function, xdata, ydata, p0=params_0, ftol=ftol, maxfev=maxfev)

    return params_opt

def get_threshold(df):
    """Returns a threshold estimate from fit."""
    
     # for L in np.unique(df['L'] == L):
     #        df = df.loc[df['L']==L]
    xdata = df[['L', 'p']].T.values
    ydata = df[['ler']].values.flatten()
    fit_data = get_fit_params(xdata, ydata, params_0=[1,1,1,1,1])
        
    return fit_data
    
def fit_function(x_data, *params):
    """ Fitting function. """
    d, p = x_data
    p_th, nu, A, B, C = params
    x = (p - p_th)*d**nu
    
    return A * x**2 + B*x + C

def get_fit_data(d, p_list, params):
    """Fitting function."""

    p_th, nu, A, B, C = params
    return [A*((p - p_th)*d**nu)**2 + B*((p - p_th)*d**nu) + C for p in p_list]

    
# Plotting tools

def plot_threshold(df, log=False):
    
    color = cm.rainbow(np.linspace(0, 1, len(np.unique(df['d']))))

    sizes = np.unique(df['d'])
    xdata = df[['d', 'p_pauli']].T.values
    ydata = df[['ler']].values.flatten()
    fit_data = get_fit_params(xdata, ydata, params_0=[1,1,1,1,1])

    for d, c in zip(sizes, color):

        dfd = df.loc[df['d'] == d]
        plt.plot(dfd['p_pauli'], dfd['ler'],'.', c=c, label=f'd = {d}')
        plt.plot(dfd['p_pauli'], get_fit_data(d, dfd['p_pauli'], fit_data), '--', c=c)
    if log == True:
        plt.loglog()
    plt.axvline(x=fit_data[0], linewidth=1, color='k', label=f'Th. fit :{round(fit_data[0]*100, 2)} %')
    plt.grid(color='gray', linestyle='dashed')
    plt.xlabel('per')
    plt.ylabel('ler')
    # plt.title(r'Logical error rate ')
    plt.savefig('data/noisy_all_p_subsystem_surface_code.png')
    plt.legend()