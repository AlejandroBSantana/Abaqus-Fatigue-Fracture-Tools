# -*- coding: utf-8 -*-
"""
Abaqus Python Script for Multiaxial Fatigue Analysis (Goodman Criterion). Abaqus version 6.14.
Based on the Maximum Principal Stress / Signed von Mises equivalent stress approach.
Developed for the Doctoral Thesis on Hyper-flexible Thoracic Implants.
"""

from odbAccess import *
import numpy as np
import time
import csv

# =====================================================================
# 1. MATHEMATICAL FUNCTIONS (SIGNED VON MISES & GOODMAN CRITERION)
# =====================================================================

def calculate_signed_mises(abaqus_tensor_list):
    """
    Computes the signed von Mises equivalent stress for a given stress history.
    The sign is determined by the hydrostatic pressure (positive for tension, 
    negative for compression).
    
    Args:
        abaqus_tensor_list (list): List of stress tensors (S11, S22, S33, S12, S13, S23) 
                                   over the load increments.
                                   
    Returns:
        np.array: Array containing the signed von Mises stress for each time increment.
    """
    num_increments = len(abaqus_tensor_list)
    signed_mises_history = np.zeros(num_increments)
    
    for t, tensor in enumerate(abaqus_tensor_list):
        s11, s22, s33, s12, s13, s23 = tensor
        
        # Calculate standard von Mises equivalent stress
        mises = np.sqrt(0.5 * ((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2 + 
                               6.0 * (s12**2 + s13**2 + s23**2)))
        
        # Calculate hydrostatic pressure to determine the stress state sign
        hydrostatic_pressure = (s11 + s22 + s33) / 3.0
        sign = 1.0 if hydrostatic_pressure >= 0 else -1.0
        
        signed_mises_history[t] = sign * mises
        
    return signed_mises_history

def evaluate_goodman_criterion(signed_mises_history, ultimate_strength, fatigue_limit):
    """
    Evaluates the Modified Goodman fatigue criterion for a given stress history.
    
    Args:
        signed_mises_history (np.array): Signed equivalent stress history over one load cycle.
        ultimate_strength (float): Ultimate tensile strength of the material (S_ut).
        fatigue_limit (float): Endurance limit of the material for fully reversed loading (S_e).
        
    Returns:
        tuple: (Goodman fatigue index, stress amplitude, mean stress, max stress, min stress)
    """
    sigma_max = np.max(signed_mises_history)
    sigma_min = np.min(signed_mises_history)
    
    stress_amplitude = (sigma_max - sigma_min) / 2.0
    mean_stress = (sigma_max + sigma_min) / 2.0
    
    # Modified Goodman criterion (ignores the infinite benefit of compressive mean stress)
    if mean_stress <= 0:
        fatigue_index = stress_amplitude / fatigue_limit
    else:
        fatigue_index = (stress_amplitude / fatigue_limit) + (mean_stress / ultimate_strength)
        
    return fatigue_index, stress_amplitude, mean_stress, sigma_max, sigma_min

# =====================================================================
# 2. ABAQUS MAIN ROUTINE & CSV EXPORT
# =====================================================================

def execute_goodman_analysis():
    """
    Main routine that opens the Abaqus ODB file, extracts the stress history at the 
    integration points of a specific element set, evaluates the Goodman criterion, 
    writes the fatigue index back into the ODB, and exports the most critical data to a CSV.
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
    
    print("Opening ODB file for Modified Goodman fatigue analysis...")
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
    print("Extracting stress history from integration points...")
    
    for frame in step.frames:
        # Ignore custom frames generated by previous fatigue analyses (e.g., Dang Van)
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
    
    print("Evaluating Modified Goodman Criterion...")
    
    for label, ips_data in stress_history.items():
        worst_fatigue_index = 0.0
        worst_ip_data = None
        
        for ip, temporal_tensors in ips_data.items():
            signed_mises = calculate_signed_mises(temporal_tensors)
            fatigue_index, s_a, s_m, s_max, s_min = evaluate_goodman_criterion(signed_mises, S_UT, S_E)
            
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

    print("\nWriting fatigue contour fields to the ODB (Frame 902)...")
    new_frame = step.Frame(incrementNumber=902, frameValue=2.0, description='Fatigue - Modified Goodman')
    fatigue_field = new_frame.FieldOutput(name='GOODMAN_INDEX', description='Goodman Damage Index', type=SCALAR)
    fatigue_field.addData(position=CENTROID, instance=instance, labels=element_labels, data=fatigue_results)
    
    odb.save()
    odb.close()
    
    print("\nGenerating CSV report (Results_Goodman.csv)...")
    with open('Results_Goodman.csv', mode='w') as file:
        writer = csv.writer(file, delimiter=';')
        
        # Double header for plotting two series in data analysis software (e.g., Excel/MATLAB)
        writer.writerow(['--- GOODMAN LIMIT LINE ---', '', '', '--- TOP 10 CRITICAL ELEMENTS ---'])
        writer.writerow(['Mean_Stress_Sm_MPa', 'Limit_Amplitude_Sa_MPa', '', 'Element_Label', 'Mean_Stress_Sm_MPa', 'Stress_Amplitude_Sa_MPa', 'Fatigue_Index_D'])
        
        # Generate 11 points to plot the theoretical Goodman limit line (from 0 to S_UT)
        sm_step = S_UT / 10.0
        for i in range(11):
            line_sm = i * sm_step
            line_sa = S_E * (1.0 - (line_sm / S_UT))
            
            # Print limit line coordinates alongside the critical elements data
            if i < len(top_critical_elements):
                el = top_critical_elements[i]
                writer.writerow([round(line_sm, 2), round(line_sa, 2), '', el['label'], round(el['Sm'], 4), round(el['Sa'], 4), round(el['Damage_Index'], 4)])
            else:
                writer.writerow([round(line_sm, 2), round(line_sa, 2), '', '', '', '', ''])

    end_time = time.time()
    print("\nAnalysis completed successfully. Fatigue contours are available in the ODB.")
    print("Total execution time: {0:.2f} seconds.".format(end_time - start_time))

if __name__ == '__main__':
    execute_goodman_analysis()