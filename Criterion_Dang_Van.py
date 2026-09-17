# -*- coding: utf-8 -*-
"""
Abaqus Python Script for Multiaxial Fatigue Analysis (Dang Van Criterion). Abaqus version 6.14.
Based on a purely NumPy-driven Badoiu-Clarkson geometric algorithm.
Compatible with Abaqus Python 2.7 environment.
Developed for the Doctoral Thesis on Hyper-flexible Thoracic Implants.
"""

from odbAccess import *
import numpy as np
import time
import csv

# =====================================================================
# 1. PURE NUMPY MATHEMATICAL FUNCTIONS (ABAQUS 6.14 / PYTHON 2.7 COMPATIBLE)
# =====================================================================

def extract_hydrostatic_and_deviatoric_tensors(abaqus_tensor_list):
    """
    Decomposes the Abaqus stress tensors into their hydrostatic pressure 
    and deviatoric stress tensor components for each time increment.
    
    Args:
        abaqus_tensor_list (list): Stress tensors (S11, S22, S33, S12, S13, S23).
                                   
    Returns:
        tuple: (hydrostatic_pressure_history, deviatoric_tensor_history)
    """
    num_increments = len(abaqus_tensor_list)
    hydrostatic_history = np.zeros(num_increments)
    deviatoric_history = np.zeros((num_increments, 3, 3))
    identity_matrix = np.eye(3)
    
    for t, tensor_abq in enumerate(abaqus_tensor_list):
        s11, s22, s33, s12, s13, s23 = tensor_abq
        stress_tensor = np.array([
            [s11, s12, s13],
            [s12, s22, s23],
            [s13, s23, s33]
        ])
        
        hydrostatic_pressure = np.trace(stress_tensor) / 3.0
        hydrostatic_history[t] = hydrostatic_pressure
        deviatoric_history[t] = stress_tensor - (hydrostatic_pressure * identity_matrix)
        
    return hydrostatic_history, deviatoric_history

def tensor_to_mandel(stress_tensor):
    """
    Converts a 3x3 symmetric tensor into a 6D Mandel notation vector.
    """
    return np.array([
        stress_tensor[0, 0], 
        stress_tensor[1, 1], 
        stress_tensor[2, 2], 
        np.sqrt(2) * stress_tensor[0, 1], 
        np.sqrt(2) * stress_tensor[0, 2], 
        np.sqrt(2) * stress_tensor[1, 2]
    ])

def mandel_to_tensor(mandel_vector):
    """
    Converts a 6D Mandel notation vector back into a 3x3 symmetric tensor.
    """
    return np.array([
        [mandel_vector[0],              mandel_vector[3]/np.sqrt(2), mandel_vector[4]/np.sqrt(2)],
        [mandel_vector[3]/np.sqrt(2),   mandel_vector[1],            mandel_vector[5]/np.sqrt(2)],
        [mandel_vector[4]/np.sqrt(2),   mandel_vector[5]/np.sqrt(2), mandel_vector[2]           ]
    ])

def calculate_residual_stress_tensor_badoiu(deviatoric_history, iterations=1000):
    """
    Executes the Badoiu-Clarkson geometric algorithm to find the microscopic 
    residual stress tensor (rho_star) that stabilizes the elastic shakedown state.
    This implementation avoids SciPy dependencies to maintain Abaqus compatibility.
    
    Args:
        deviatoric_history (np.array): History of deviatoric stress tensors.
        iterations (int): Number of iterations for convergence.
        
    Returns:
        np.array: 3x3 microscopic residual stress tensor (rho_star).
    """
    mandel_history = np.array([tensor_to_mandel(S) for S in deviatoric_history])
    center = np.mean(mandel_history, axis=0)
    
    for i in range(1, iterations + 1):
        squared_distances = np.sum((mandel_history - center)**2, axis=1)
        farthest_idx = np.argmax(squared_distances)
        farthest_point = mandel_history[farthest_idx]
        
        step_size = 1.0 / (i + 1)
        center = center + step_size * (farthest_point - center)
        
    rho_star = mandel_to_tensor(-center)
    return rho_star

def evaluate_dang_van_criterion(deviatoric_history, hydrostatic_history, rho_star, const_a, const_b):
    """
    Evaluates the multiaxial Dang Van fatigue criterion for the load cycle.
    
    Args:
        deviatoric_history (np.array): Deviatoric stress tensor history.
        hydrostatic_history (np.array): Hydrostatic pressure history.
        rho_star (np.array): Microscopic residual stress tensor.
        const_a (float): Material-specific Dang Van constant A.
        const_b (float): Material-specific Dang Van constant B (MPa).
        
    Returns:
        tuple: (Maximum fatigue index, microscopic shear stress history)
    """
    max_damage_index = 0.0
    shear_stress_history = [] 
    
    for t in range(len(hydrostatic_history)):
        local_stress = deviatoric_history[t] + rho_star
        eigenvalues = np.linalg.eigvalsh(local_stress)
        principal_stresses = np.sort(eigenvalues)[::-1]
        
        micro_shear_stress = (principal_stresses[0] - principal_stresses[2]) / 2.0
        shear_stress_history.append(micro_shear_stress)
        
        current_damage_index = (micro_shear_stress + const_a * hydrostatic_history[t]) / const_b
        if current_damage_index > max_damage_index:
            max_damage_index = current_damage_index
            
    return max_damage_index, shear_stress_history

# =====================================================================
# 2. ABAQUS MAIN ROUTINE (DATA EXTRACTION AND ODB WRITING)
# =====================================================================

def execute_dang_van_analysis():
    """
    Main execution routine for the Dang Van multiaxial fatigue evaluation.
    Extracts stress data from the ODB, calculates the fatigue index for all 
    elements in a specific set, writes the contours to the ODB, and exports 
    the critical element path to a CSV for scientific publication plotting.
    """
    # --- MODEL CONFIGURATION ---
    ODB_PATH = 'breathing.odb' 
    INSTANCE_NAME = 'PROTESIS3_PERFORADA'
    STEP_NAME = 'Desplazamientos'         
    CRITICAL_SET_NAME = 'REFINAMIENTO' 
    
    # --- MATERIAL PROPERTIES (Ti-6Al-4V ELI LPBF) ---
    CONSTANT_A = 0.24 
    CONSTANT_B = 195.0 # MPa
    # ----------------------------------------------
    
    print("Opening ODB file for Dang Van multiaxial fatigue analysis...")
    start_time = time.time()
    
    try:
        odb = openOdb(path=ODB_PATH, readOnly=False)
    except Exception as e:
        print("[ERROR] Could not open the ODB file. Please close it in Abaqus/Viewer if open.")
        return

    try:
        instance = odb.rootAssembly.instances[INSTANCE_NAME]
        step = odb.steps[STEP_NAME]
    except KeyError:
        print("[ERROR] Instance or Step name not found in the ODB.")
        odb.close()
        return

    if CRITICAL_SET_NAME in instance.elementSets.keys():
        element_set = instance.elementSets[CRITICAL_SET_NAME]
        print(" -> Set '" + CRITICAL_SET_NAME + "' found in Instance.")
    elif CRITICAL_SET_NAME in odb.rootAssembly.elementSets.keys():
        element_set = odb.rootAssembly.elementSets[CRITICAL_SET_NAME]
        print(" -> Set '" + CRITICAL_SET_NAME + "' found in Assembly.")
    else:
        print("\n[FATAL ERROR] Element Set not found: " + CRITICAL_SET_NAME)
        print("Available Instance Sets: " + str(instance.elementSets.keys()))
        print("Available Assembly Sets: " + str(odb.rootAssembly.elementSets.keys()))
        odb.close()
        return
    
    stress_history = {}
    
    print("Extracting stress history (S) across all increments...")
    for frame in step.frames:
        if frame.incrementNumber > 900: 
            continue
            
        stress_field = frame.fieldOutputs['S']
        stress_values = stress_field.getSubset(region=element_set).values
        
        for value in stress_values:
            label = value.elementLabel
            ip = value.integrationPoint
            
            if label not in stress_history:
                stress_history[label] = {}
            if ip not in stress_history[label]:
                stress_history[label][ip] = []
                
            stress_history[label][ip].append(value.data)

    element_labels = []
    fatigue_results = []
    
    global_worst_element = None
    global_worst_damage_index = -1.0
    critical_export_data = {}
    
    num_elements = len(stress_history)
    print("Evaluating Dang Van Criterion for {0} C3D10 elements...".format(num_elements))
    
    counter = 0
    for label, ips_data in stress_history.items():
        worst_element_dv_index = 0.0
        worst_element_ip_data = None
        
        for ip, temporal_tensors in ips_data.items():
            hydro_hist, dev_hist = extract_hydrostatic_and_deviatoric_tensors(temporal_tensors)
            rho_star = calculate_residual_stress_tensor_badoiu(dev_hist)
            dv_index, tau_hist = evaluate_dang_van_criterion(dev_hist, hydro_hist, rho_star, CONSTANT_A, CONSTANT_B)
            
            if dv_index > worst_element_dv_index:
                worst_element_dv_index = dv_index
                worst_element_ip_data = {'ph': hydro_hist, 'tau': tau_hist}
                
        element_labels.append(label)
        fatigue_results.append((worst_element_dv_index, ))
        
        if worst_element_dv_index > global_worst_damage_index:
            global_worst_damage_index = worst_element_dv_index
            global_worst_element = label
            critical_export_data = worst_element_ip_data
        
        counter += 1
        if counter % 50 == 0:
            print("  Progress: {0}/{1} elements processed".format(counter, num_elements))

    print("\nWriting fatigue contour fields back to the ODB (Centroid)...")
    new_frame = step.Frame(incrementNumber=903, frameValue=3.0, description='Fatigue - Dang Van')
    fatigue_field = new_frame.FieldOutput(name='DANG_VAN_INDEX', description='Dang Van Damage Index', type=SCALAR)
    fatigue_field.addData(position=CENTROID, instance=instance, labels=element_labels, data=fatigue_results)
    
    odb.save()
    odb.close()
    
    # --- CRITICAL ELEMENT EXPORT FOR PUBLICATION PLOTS ---
    print("\nExporting critical element load path (Label: {0}) for publication plots...".format(global_worst_element))
    csv_filename = 'Critical_Element_DangVan_Path.csv'
    with open(csv_filename, mode='wb') as csv_file:
        writer = csv.writer(csv_file, delimiter=';')
        writer.writerow(['Increment', 'Hydrostatic_Pressure_Ph_MPa', 'Microscopic_Shear_Stress_Tau_MPa', 'Material_Limit_Line'])
        
        for t in range(len(critical_export_data['ph'])):
            ph_val = critical_export_data['ph'][t]
            tau_val = critical_export_data['tau'][t]
            limit_val = CONSTANT_B - (CONSTANT_A * ph_val) 
            writer.writerow([t+1, round(ph_val, 4), round(tau_val, 4), round(limit_val, 4)])
            
    print("File {0} successfully generated.".format(csv_filename))
    # -------------------------------------------------------

    end_time = time.time()
    print("\n=== ANALYSIS SUCCESSFULLY COMPLETED IN {0} SECONDS ===".format(round(end_time - start_time, 2)))
    print("Open the ODB in Abaqus Viewer to visualize the DANG_VAN_INDEX contour.")

if __name__ == '__main__':
    execute_dang_van_analysis()