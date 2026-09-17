# -*- coding: utf-8 -*-
"""
Abaqus Python Script for Multiaxial Fatigue Analysis (Goodman Criterion). Abaqus version 6.14.
Based on the Maximum Principal Stress (Sigma 1) approach for Mode I crack opening.
Developed for the Doctoral Thesis on Hyper-flexible Thoracic Implants.
"""

from odbAccess import *
import numpy as np
import time
import csv

# =====================================================================
# 1. MATHEMATICAL FUNCTIONS (MAXIMUM PRINCIPAL STRESS & GOODMAN)
# =====================================================================

def calculate_sigma1_history(abaqus_tensor_list):
    """
    Extracts the time history of the maximum principal stress (sigma_1) 
    for a given integration point. This stress is primarily responsible 
    for Mode I crack opening.
    
    Args:
        abaqus_tensor_list (list): List of stress tensors (S11, S22, S33, S12, S13, S23).
                                   
    Returns:
        np.array: Array containing the maximum principal stress for each time increment.
    """
    history = np.zeros(len(abaqus_tensor_list))
    
    for i, tensor_abq in enumerate(abaqus_tensor_list):
        s11, s22, s33, s12, s13, s23 = tensor_abq
        
        stress_tensor = np.array([
            [s11, s12, s13],
            [s12, s22, s23],
            [s13, s23, s33]
        ])
        
        # Calculate eigenvalues of the symmetric stress tensor
        principal_stresses = np.linalg.eigvalsh(stress_tensor)
        
        # np.linalg.eigvalsh returns sorted eigenvalues: sigma3 <= sigma2 <= sigma1
        history[i] = principal_stresses[2]
        
    return history

def evaluate_goodman_criterion(sigma1_history, ultimate_strength, fatigue_limit):
    """
    Evaluates the Modified Goodman fatigue criterion using the maximum principal stress history.
    
    Args:
        sigma1_history (np.array): Sigma_1 stress history over one load cycle.
        ultimate_strength (float): Ultimate tensile strength of the material (S_ut).
        fatigue_limit (float): Endurance limit of the material for fully reversed loading (S_e).
        
    Returns:
        tuple: (Goodman fatigue index, stress amplitude, mean stress, max stress, min stress)
    """
    sigma_max = np.max(sigma1_history)
    sigma_min = np.min(sigma1_history)
    
    stress_amplitude = (sigma_max - sigma_min) / 2.0
    mean_stress = (sigma_max + sigma_min) / 2.0
    
    # Modified Goodman criterion (conservative assumption for compressive mean stress)
    if mean_stress < 0.0:
        fatigue_index = stress_amplitude / fatigue_limit
    else:
        fatigue_index = (stress_amplitude / fatigue_limit) + (mean_stress / ultimate_strength)
        
    return fatigue_index, stress_amplitude, mean_stress, sigma_max, sigma_min

# =====================================================================
# 2. ABAQUS MAIN ROUTINE & CSV EXPORT
# =====================================================================

def execute_goodman_tpmax_analysis():
    """
    Main routine that opens the Abaqus ODB file, extracts the stress history at the 
    integration points of a specific element set, evaluates the Goodman criterion 
    based on Maximum Principal Stress, writes the fatigue index back into the ODB, 
    and exports the most critical data to a CSV file.
    """
    # --- MODEL CONFIGURATION ---
    ODB_PATH = 'breathing.odb' 
    INSTANCE_NAME = 'PROTESIS3_PERFORADA' 
    STEP_NAME = 'Desplazamientos'         
    CRITICAL_SET_NAME = 'REFINAMIENTO'
    
    # --- MATERIAL PROPERTIES (Ti-6Al-4V ELI LPBF) ---
    S_UT = 1150.88  # Ultimate Tensile Strength (MPa)
    S_E = 336.0     # Fatigue Limit at R=-1 (MPa)
    # ----------------------------------------------
    
    print("Opening ODB file for Goodman Analysis (Maximum Principal Stress)...")
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
    
    # Locate the critical element set
    if CRITICAL_SET_NAME in instance.elementSets.keys():
        element_set = instance.elementSets[CRITICAL_SET_NAME]
    elif CRITICAL_SET_NAME in odb.rootAssembly.elementSets.keys():
        element_set = odb.rootAssembly.elementSets[CRITICAL_SET_NAME]
    else:
        print("[ERROR] Element set not found: " + CRITICAL_SET_NAME)
        odb.close()
        return
    
    stress_history = {}
    print("Extracting stress tensors from the model integration points...")
    
    for frame in step.frames:
        # Ignore custom frames generated by previous fatigue analyses
        if frame.incrementNumber > 900: 
            continue
            
        stress_field = frame.fieldOutputs['S']
        stress_values = stress_field.getSubset(region=element_set).values
        
        for value in stress_values:
            element_label = value.elementLabel
            integration_point = value.integrationPoint
            
            if element_label not in stress_history:
                stress_history[element_label] = {}
            if integration_point not in stress_history[element_label]:
                stress_history[element_label][integration_point] = []
                
            stress_history[element_label][integration_point].append(value.data)

    element_labels = []
    fatigue_results = []
    element_data_list = []
    
    print("Evaluating Maximum Principal Stress Goodman Criterion...")
    
    for label, ips_data in stress_history.items():
        worst_fatigue_index = 0.0
        worst_ip_data = None
        
        for ip, temporal_tensors in ips_data.items():
            sigma1_history = calculate_sigma1_history(temporal_tensors)
            fatigue_index, s_a, s_m, s_max, s_min = evaluate_goodman_criterion(sigma1_history, S_UT, S_E)
            
            # Store the worst integration point for the current element
            if fatigue_index > worst_fatigue_index:
                worst_fatigue_index = fatigue_index
                worst_ip_data = {'label': label, 'Sm': s_m, 'Sa': s_a, 'Damage_Index': fatigue_index}
                
        element_labels.append(label)
        fatigue_results.append((worst_fatigue_index, ))
        element_data_list.append(worst_ip_data)

    # Sort elements by descending fatigue damage index
    element_data_list.sort(key=lambda x: x['Damage_Index'], reverse=True)
    top_critical_elements = element_data_list[:10]

    print("\nWriting fatigue contour fields to the ODB (Frame 901)...")
    new_frame = step.Frame(incrementNumber=901, frameValue=1.0, description='Fatigue - Goodman TPMAX')
    fatigue_field = new_frame.FieldOutput(name='GOODMAN_INDEX', description='Goodman Damage Index', type=SCALAR)
    fatigue_field.addData(position=CENTROID, instance=instance, labels=element_labels, data=fatigue_results)
    
    odb.save()
    odb.close()
    
    print("\nGenerating CSV report (Results_Goodman_TPMAX.csv)...")
    with open('Results_Goodman_TPMAX.csv', mode='w') as file:
        writer = csv.writer(file, delimiter=';')
        
        writer.writerow(['--- GOODMAN LIMIT LINE ---', '', '', '--- TOP 10 CRITICAL ELEMENTS ---'])
        writer.writerow(['Mean_Stress_Sm_MPa', 'Limit_Amplitude_Sa_MPa', '', 'Element_Label', 'Mean_Stress_Sm_MPa', 'Stress_Amplitude_Sa_MPa', 'Fatigue_Index_D'])
        
        sm_step = S_UT / 10.0
        for i in range(11):
            line_sm = i * sm_step
            line_sa = S_E * (1.0 - (line_sm / S_UT))
            
            if i < len(top_critical_elements):
                el = top_critical_elements[i]
                writer.writerow([round(line_sm, 2), round(line_sa, 2), '', el['label'], round(el['Sm'], 4), round(el['Sa'], 4), round(el['Damage_Index'], 4)])
            else:
                writer.writerow([round(line_sm, 2), round(line_sa, 2), '', '', '', '', ''])

    end_time = time.time()
    print("\nAnalysis completed successfully in {0:.2f} seconds.".format(end_time - start_time))

if __name__ == '__main__':
    execute_goodman_tpmax_analysis()