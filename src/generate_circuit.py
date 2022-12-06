"""
This file contains a function to generate the rotated subsystem surface code for a memory experiment.
"""

from typing import Set, List, Dict
from dataclasses import dataclass
from CircuitGeneratorParams import *

@dataclass
class SubsystemSurfaceCode:
    distance : int 
    
    def get_subsystem_surface_code_coordinates(self, distance) -> List[complex]:
        """Return a list of relevant coordinates in the surface code."""

        # Place data qubits
        data_coords: Set[complex] = set()
        x_observable: List[complex] = []
        z_observable: List[complex] = []
        # Place measurement qubits.
        x_ancilla_coords: Set[complex] = set()
        z_ancilla_coords: Set[complex] = set()

        for x in range(2*distance+1):
            for y in range(2*distance+1):
                if x % 4 == 0 and y % 4 == 2 or x % 4 == 2 and y % 4 == 0 :
                    if 0 <= y <= 2*distance+1 and 0<x<2*distance:
                        x_ancilla_coords.add(x+1j*y)
                if x % 4 == 0 and y % 4 == 0 or x % 4 == 2 and y % 4 == 2 :
                    if 0 <= x <=2*distance+1 and 0<y<2*distance:
                        z_ancilla_coords.add(x+1j*y)
                if x % 4 == 1 and y % 4 == 1:
                    data_coords.add(x+1j*y)
                if x % 4 == 3 and y % 4 == 1 or y % 4 == 3 and x % 4 == 1:
                    data_coords.add(x+1j*y)
        for c in data_coords:
            if c.real == 1:
                x_observable.append(c)
            elif c.imag == 1:
                z_observable.append(c)

        return data_coords, x_ancilla_coords, z_ancilla_coords, x_observable, z_observable 

    @property 
    def data_coords(self) -> List[complex]:
        data_coords, x_ancilla_coords, z_ancilla_coords, x_observable, z_observable = self.get_subsystem_surface_code_coordinates(self.distance)
        return data_coords
    
    @property 
    def x_ancilla_coords(self) -> List[complex]:
        data_coords, x_ancilla_coords, z_ancilla_coords, x_observable, z_observable = self.get_subsystem_surface_code_coordinates(self.distance)
        return x_ancilla_coords
        
    @property 
    def z_ancilla_coords(self) -> List[complex]:
        data_coords, x_ancilla_coords, z_ancilla_coords, x_observable, z_observable = self.get_subsystem_surface_code_coordinates(self.distance)
        return z_ancilla_coords
    
    @dataclass 
    class TriangleType:
        pauli: str
        ancilla: complex
        round_order: int
        horizontal_link : complex
        vertical_link : complex
        diagonal_link : complex
        link_order : Dict[int, complex]
        gauge_pair : complex 
        
    @property 
    def TriangleDict(self):
        
        triangleDict = {}

        for q in self.x_ancilla_coords:
            if q.real % 4 == 2 and q.imag % 4 == 0:
                triangleDict[q] = self.TriangleType(pauli="X" , ancilla = q, round_order=0, horizontal_link=-1-1j, vertical_link=1+1j, diagonal_link=-1+1j,
                                                    link_order={0: None, 1: -1+1j, 2: -1-1j, 3: 1+1j}, gauge_pair=q+2-2j)
               
            elif q.real % 4 == 0 and q.imag % 4 == 2:
                triangleDict[q] = self.TriangleType(pauli="X", ancilla = q, round_order=3, horizontal_link=1+1j, vertical_link=-1-1j, diagonal_link=1-1j,
                                                    link_order={0: 1-1j, 1: 1+1j, 2: -1-1j, 3: None}, gauge_pair=None)
        
        for q in self.z_ancilla_coords:
            if q.real % 4 == 2 and q.imag % 4 == 2:
                triangleDict[q] = self.TriangleType(pauli="Z", ancilla = q, round_order=2, horizontal_link=-1+1j, vertical_link=1-1j, diagonal_link=-1-1j,
                                                    link_order={0: -1+1j, 1: 1-1j, 2: None, 3: -1-1j}, gauge_pair=q+2+2j)
                
            elif q.real % 4 == 0 and q.imag % 4 == 0:
                triangleDict[q] = self.TriangleType(pauli="Z", ancilla = q, round_order=1, horizontal_link=1-1j, vertical_link=-1+1j, diagonal_link=1+1j,
                                                    link_order={0: -1+1j, 1: None, 2: 1+1j, 3: 1-1j}, gauge_pair=None)

        return triangleDict 
    
    @property
    def measure_at_round_coords(self):
        
        m : List[List[complex]] = [[], [], [], []]
        for measure in self.TriangleDict:
            k = self.TriangleDict[measure].round_order
            m[k].append(measure)
        return m
    
    
    @property
    def m2shape(self):
        
        m2shape: Dict[complex, str] = {}
        for m in list(self.x_ancilla_coords)+list(self.z_ancilla_coords):
            if m.real in [0, 2*self.distance] or m.imag in [0, 2*self.distance]:
                m2shape[m] = 'stabilizer'
            else:
                m2shape[m] = 'gauge'
        return m2shape
            

def generate_rotated_subsystem_surface_code(params: CircuitGenParameters, codeparams: SubsystemSurfaceCode, experiment_basis: str) -> stim.Circuit:
    """Generate the circuit for the subsystem surface code as per Bravyi's paper."""
    
    if params.rounds < 1:
        raise ValueError("Need rounds >= 1")
    if params.distance < 3:
        raise ValueError("Need a distance >= 2")
    
    # Create qubit coords
    data_coords, x_ancilla_coords, z_ancilla_coords, x_observable, z_observable = codeparams.get_subsystem_surface_code_coordinates(params.distance)
    
    # Index the measurement qubits and data qubits.
    p2q: Dict[complex, int] = {}
    for q in data_coords:
        p2q[q] = coord_to_idx(q, params.distance)

    for q in x_ancilla_coords:
        p2q[q] = coord_to_idx(q, params.distance)

    for q in z_ancilla_coords:
        p2q[q] = coord_to_idx(q, params.distance)

    q2p: Dict[int, complex] = {v: k for k, v in p2q.items()}
        
    data_qubits = [p2q[q] for q in data_coords]
    measurement_qubits = [p2q[q] for q in x_ancilla_coords]
    measurement_qubits += [p2q[q] for q in z_ancilla_coords]
    x_ancilla_qubits = [p2q[q] for q in x_ancilla_coords]
    z_ancilla_qubits = [p2q[q] for q in z_ancilla_coords]

    all_qubits: List[int] = []
    all_qubits += data_qubits + measurement_qubits

    all_qubits.sort()
    data_qubits.sort()
    measurement_qubits.sort()
    x_ancilla_qubits.sort()
    z_ancilla_qubits.sort()
        
    # Reverse index the measurement order used for defining detectors
    data_coord_to_order: Dict[complex, int] = {}
    measure_coord_to_order: Dict[complex, int] = {}
        
    for q in data_qubits:
        data_coord_to_order[q2p[q]] = len(data_coord_to_order)
        
    for q in measurement_qubits:
        measure_coord_to_order[q2p[q]] = len(measure_coord_to_order)
    
    # Indicate which ancillas are of same type as the experiment basis 
    chosen_basis_measure_coords = x_ancilla_coords if experiment_basis == "X" else z_ancilla_coords
    chosen_basis_observable = x_observable if experiment_basis == "X" else z_observable 

    
    # List out CNOT gate targets in the bulk 
    #---------------------------------------

    cnot_targets: List[List[int]] = [[], [], [], []]
    
    for k in range(4):
                
        for measure in sorted(x_ancilla_coords, key=lambda c: (c.real, c.imag)):
            tri_op = codeparams.TriangleDict[measure]

            if tri_op.round_order is not k:

                data = measure + tri_op.link_order[k]
                if data in p2q:
                    # print(p2q[measure], p2q[data])
                    cnot_targets[k].append(p2q[measure])
                    cnot_targets[k].append(p2q[data])

        for measure in sorted(z_ancilla_coords, key=lambda c: (c.real, c.imag)):
            
            tri_op = codeparams.TriangleDict[measure]
                
            if tri_op.round_order is not k:

                data = measure + tri_op.link_order[k]
                if data in p2q:
                    cnot_targets[k].append(p2q[data])
                    cnot_targets[k].append(p2q[measure])   

    # List out CNOT gate targets at start 
    #------------------------------------
    start_cnot_targets: List[List[int]] = [[], [], [], []]
    
    for k in range(2):
        
        for measure in sorted(x_ancilla_coords, key=lambda c: (c.real, c.imag)):
            tri_op = codeparams.TriangleDict[measure]

            if tri_op.round_order % 3 == 0 and tri_op.round_order != k:

                data = measure + tri_op.link_order[k]
                if data in p2q:
                    start_cnot_targets[k].append(p2q[measure])
                    start_cnot_targets[k].append(p2q[data])

        for measure in sorted(z_ancilla_coords, key=lambda c: (c.real, c.imag)):
            tri_op = codeparams.TriangleDict[measure]

            if tri_op.round_order % 3 == 0 and tri_op.round_order != k:

                data = measure + tri_op.link_order[k]
                if data in p2q:
                    start_cnot_targets[k].append(p2q[data])
                    start_cnot_targets[k].append(p2q[measure])  
                    
    start_cnot_targets[2] = cnot_targets[2].copy()
    start_cnot_targets[3] = cnot_targets[3].copy()
    
    # List out CNOT gate targets at end 
    #----------------------------------
    end_cnot_targets: List[List[int]] = [[], [], [], []]
    end_cnot_targets[0] = cnot_targets[0].copy()
    
    for k in range(1,4):
        
        for measure in sorted(x_ancilla_coords, key=lambda c: (c.real, c.imag)):
            tri_op = codeparams.TriangleDict[measure]
            if k < tri_op.round_order:
                data = measure + tri_op.link_order[k]
                if data in p2q:
                    end_cnot_targets[k].append(p2q[measure])
                    end_cnot_targets[k].append(p2q[data])
                    
        for measure in sorted(z_ancilla_coords, key=lambda c: (c.real, c.imag)):
            tri_op = codeparams.TriangleDict[measure]
            if k < tri_op.round_order:
                data = measure + tri_op.link_order[k]
                if data in p2q:
                    end_cnot_targets[k].append(p2q[data])
                    end_cnot_targets[k].append(p2q[measure])      
    
    # Bulk CNOT and measurement cycle 
    #--------------------------------
    measurement_series : List[int] = []

    bulk_cycle_actions = stim.Circuit()
    params.append_begin_round_tick(bulk_cycle_actions, data_qubits)  
    
    for step, targets in enumerate(cnot_targets):
        
        params.append_unitary_2(bulk_cycle_actions, "CNOT", targets)
        if step in [0,3]: 
            params.append_measure_reset(bulk_cycle_actions, [p2q[i] for i in codeparams.measure_at_round_coords[step]], basis='X')
            measurement_series += [p2q[i] for i in codeparams.measure_at_round_coords[step]]
        elif step in [1,2]: 
            params.append_measure_reset(bulk_cycle_actions, [p2q[i] for i in codeparams.measure_at_round_coords[step]])
            measurement_series += [p2q[i] for i in codeparams.measure_at_round_coords[step]]
        bulk_cycle_actions.append_operation("TICK", [])

    # -------- BUILD CIRCUIT --------

    subsystem_measure_coord_to_order: Dict[complex, int] = {}

    for q in measurement_series:
        subsystem_measure_coord_to_order[q2p[q]] = len(subsystem_measure_coord_to_order)
    
    # Build the starting cycle  
    #-------------------------
    head = stim.Circuit()
    if params.show_qubit_coordinates == True:
        for k, v in sorted(q2p.items()):
            head.append_operation("QUBIT_COORDS", [k], [v.real, v.imag])
        
    # Reset circuit 
    params.append_reset(head, data_qubits, experiment_basis)
    params.append_reset(head, z_ancilla_qubits)
    params.append_reset(head, x_ancilla_qubits, experiment_basis)
    params.append_begin_round_tick(head, data_qubits)

    for step, targets in enumerate(start_cnot_targets):
        if step == 3: 
            params.append_unitary_2(head, "CNOT", targets)
            first_measure_qubits = [p2q[i] for i in codeparams.measure_at_round_coords[step]]
            params.append_measure_reset(head, [p2q[i] for i in codeparams.measure_at_round_coords[step]], basis='X')
            head.append_operation("TICK", [])
        else:
            params.append_unitary_2(head, "CNOT", targets)
            head.append_operation("TICK", [])

    # First detectors in experiment basis
    # boundary (stabilizer detectors)
    m = len(measurement_series)
    for measure in sorted(chosen_basis_measure_coords, key=lambda c: (c.real, c.imag)):
        
        if codeparams.m2shape[measure] == 'stabilizer' and codeparams.TriangleDict[measure].round_order == 3:
            k = m - subsystem_measure_coord_to_order[measure] 
            head.append_operation(
                "DETECTOR",
                [stim.target_rec(-k)],
                [measure.real, measure.imag, 0.0]
            )
    
    # Then need to add a cycle 
    head += bulk_cycle_actions
    
    # Add remaining first experiment basis detectors 
    for measure in sorted(chosen_basis_measure_coords, key=lambda c: (c.real, c.imag)):
        if codeparams.m2shape[measure] == 'stabilizer':
            
            if codeparams.TriangleDict[measure].round_order == 3:
                head.append_operation(
                    "DETECTOR",
                    [stim.target_rec(-2*m + subsystem_measure_coord_to_order[measure]),
                    stim.target_rec(-m + subsystem_measure_coord_to_order[measure])],
                    [measure.real, measure.imag, 0.0]
                    )
            else:
                head.append_operation(
                    "DETECTOR",
                    [stim.target_rec(-m + subsystem_measure_coord_to_order[measure])],
                    [measure.real, measure.imag, 0.0]
                    )
                
        else:
            gauge_pair = codeparams.TriangleDict[measure].gauge_pair
            if gauge_pair is not None:

                head.append_operation(
                    "DETECTOR",
                    [stim.target_rec(-m + subsystem_measure_coord_to_order[measure]),
                     stim.target_rec(-2*m + subsystem_measure_coord_to_order[gauge_pair])],
                    [measure.real, measure.imag, 0.0]) 
                
    # Build the bulk 
    # --------------
    
    body = stim.Circuit()

    body += bulk_cycle_actions.copy()
    
    body.append_operation("SHIFT_COORDS", [], [0.0, 0.0, 1.0])   
    
    for m_index in measurement_series:
        m_coord = q2p[m_index]
        
        # Add stabilizer detectors (boundary) first : both X and Z 
        if codeparams.m2shape[m_coord] == 'stabilizer':

            k = m - subsystem_measure_coord_to_order[m_coord] - 1

            if not params.exclude_other_basis_detectors or m_coord in chosen_basis_measure_coords:
                body.append_operation(
                    "DETECTOR",
                    [stim.target_rec(-k - 1), stim.target_rec(-k - 1 - m)],
                    [m_coord.real, m_coord.imag, 0.0])
            
        # Add gauge detectors (boundary) first : both X and Z 
        else:
            triangle = codeparams.TriangleDict[m_coord]
            if not params.exclude_other_basis_detectors or m_coord in chosen_basis_measure_coords:
                if triangle.gauge_pair is not None:
                    if triangle.round_order == 0:
                        k1 = m - subsystem_measure_coord_to_order[m_coord] - 1
                        k2 = 2*m - subsystem_measure_coord_to_order[triangle.gauge_pair] - 1
                    else:
                        k1 = m - subsystem_measure_coord_to_order[m_coord] - 1
                        k2 = m - subsystem_measure_coord_to_order[triangle.gauge_pair] - 1
                    body.append_operation(
                        "DETECTOR",
                        [stim.target_rec(-k1 - 1), stim.target_rec(-k1 - 1 - m),
                         stim.target_rec(-k2 - 1), stim.target_rec(-k2 - 1 - m)],
                        [m_coord.real, m_coord.imag, 0.0]) 

    # Build the end  
    #---------------

    tail = stim.Circuit()
    tail.append_operation("TICK", [])

    for step, targets in enumerate(end_cnot_targets):
        
        if step == 0:
            params.append_unitary_2(tail, "CNOT", targets)
            params.append_measure_reset(tail, [p2q[i] for i in codeparams.measure_at_round_coords[step]], basis='X')
            # end_measurement_series += [p2q[i] for i in codeparams.measure_at_round_coords[step]]
            
        elif step in [1,2]: 
            params.append_unitary_2(tail, "CNOT", targets)
            params.append_measure_reset(tail, [p2q[i] for i in codeparams.measure_at_round_coords[step]])
            # end_measurement_series += [p2q[i] for i in codeparams.measure_at_round_coords[step]]
            
        else:
            params.append_measure_reset(tail, [p2q[i] for i in codeparams.measure_at_round_coords[step]], basis='X')
        tail.append_operation("TICK", [])


    # Before-final detectors
    #-----------------------

    for m_index in measurement_series:
        m_coord = q2p[m_index]
        
        # Add stabilizer detectors (boundary) first : both X and Z 
        if codeparams.m2shape[m_coord] == 'stabilizer':

            k = m - subsystem_measure_coord_to_order[m_coord] - 1

            if not params.exclude_other_basis_detectors or m_coord in chosen_basis_measure_coords:
                tail.append_operation(
                    "DETECTOR",
                    [stim.target_rec(-k - 1), stim.target_rec(-k - 1 - m)],
                    [m_coord.real, m_coord.imag, 0.0])
            
        # Add gauge detectors (boundary) first : both X and Z 
        else:
            triangle = codeparams.TriangleDict[m_coord]
            if not params.exclude_other_basis_detectors or m_coord in chosen_basis_measure_coords:
                if triangle.gauge_pair is not None:
                    if triangle.round_order == 0:
                        k1 = m - subsystem_measure_coord_to_order[m_coord] - 1
                        k2 = 2*m - subsystem_measure_coord_to_order[triangle.gauge_pair] - 1
                    else:
                        k1 = m - subsystem_measure_coord_to_order[m_coord] - 1
                        k2 = m - subsystem_measure_coord_to_order[triangle.gauge_pair] - 1
                    tail.append_operation(
                        "DETECTOR",
                        [stim.target_rec(-k1 - 1), stim.target_rec(-k1 - 1 - m),
                         stim.target_rec(-k2 - 1), stim.target_rec(-k2 - 1 - m)],
                        [m_coord.real, m_coord.imag, 0.0]) 
    # Final detectors
    #----------------

    params.append_measure(tail, data_qubits, experiment_basis)
        
    for measure in sorted(chosen_basis_measure_coords, key=lambda c: (c.real, c.imag)):
        triangle = codeparams.TriangleDict[measure]
        detectors: List[int] = []
        if codeparams.m2shape[measure] == 'stabilizer':
            for delta in [triangle.horizontal_link, triangle.vertical_link, triangle.diagonal_link]:
                data = measure + delta
                if data in data_coords:
                    detectors.append(-len(data_qubits) + data_coord_to_order[data])
            detectors.append(-len(data_qubits) - m + subsystem_measure_coord_to_order[measure])
            detectors.sort(reverse=True)
            tail.append_operation("DETECTOR", [stim.target_rec(x) for x in detectors], [measure.real, measure.imag, 1.0])
            
        else:
            triangle = codeparams.TriangleDict[measure]
            if triangle.gauge_pair is not None:
                gauge_pair = codeparams.TriangleDict[triangle.gauge_pair]
                for delta in [triangle.horizontal_link, triangle.vertical_link, triangle.diagonal_link]:
                    data = measure + delta
                    if data in data_coords:
                        detectors.append(-len(data_qubits) + data_coord_to_order[data])
                for delta in [gauge_pair.horizontal_link, gauge_pair.vertical_link, gauge_pair.diagonal_link]:  
                    data = triangle.gauge_pair + delta
                    if data in data_coords:
                        detectors.append(-len(data_qubits) + data_coord_to_order[data])
                        
                if triangle.round_order == 3:
                    detectors.append(-len(data_qubits) - 2*m + subsystem_measure_coord_to_order[measure])
                else:
                    detectors.append(-len(data_qubits) - m + subsystem_measure_coord_to_order[measure])
                    
                if gauge_pair.round_order == 3:
                    detectors.append(-len(data_qubits) - 2*m + subsystem_measure_coord_to_order[triangle.gauge_pair])
                else:
                    detectors.append(-len(data_qubits) - m + subsystem_measure_coord_to_order[triangle.gauge_pair])
                detectors.sort(reverse=True)
                tail.append_operation("DETECTOR", [stim.target_rec(x) for x in detectors], [measure.real, measure.imag, 1.0])
            
 
    # Logical observable
    # -------------------
    obs_inc: List[int] = []
    for q in chosen_basis_observable:
        obs_inc.append(-len(data_qubits) + data_coord_to_order[q])
    obs_inc.sort(reverse=True)
    tail.append_operation("OBSERVABLE_INCLUDE", [stim.target_rec(x) for x in obs_inc], 0.0)

    return head + body * (params.rounds - 1) + tail
    
# circuit = generate_rotated_subsystem_surface_code(CircuitGenParameters(distance=3, rounds=5, show_qubit_coordinates=False), SubsystemSurfaceCode(distance=3), experiment_basis='X')

# print(repr(circuit))
# sampler = circuit.compile_sampler()
# arr = sampler.sample(shots=1)
# print(repr(circuit.detector_error_model(allow_gauge_detectors=False)))