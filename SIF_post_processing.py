# -*- coding: utf-8 -*-
"""
Abaqus Python Script for Stress Intensity Factor (SIF) Post-Processing. Abaqus version 6.14.
Automatically extracts and filters fracture mechanics data (Interaction Integral contours).
Note: Currently configured for Mode I (K1). To evaluate Mode II (K2) or Mode III (K3),
simply change the target variable in the data extraction loop (e.g., replace 'K1' with 'K2').
Developed for the Doctoral Thesis on Hyper-flexible Thoracic Implants.
"""

from odbAccess import openOdb
import csv
import re

# =====================================================================
# 1. USER CONFIGURATION
# =====================================================================
ODB_NAME = 'Job-1.odb'            # Exact name of the ODB file
STEP_NAME = 'Desplazamientos'     # Name of the analysis step
CSV_FILENAME = 'Results_K1_a_010.csv' # Output file name
TOLERANCE = 0.05                  # Maximum allowed dispersion (5%)
MIN_CONTOURS = 3                  # Minimum number of grouped contours required for a valid plateau
# =====================================================================

def execute_sif_post_processing():
    print("Opening ODB and reading history output data...")
    try:
        odb = openOdb(path=ODB_NAME, readOnly=True)
    except Exception as e:
        print("[ERROR] Could not open the ODB file. Ensure it is not locked by Abaqus/Viewer.")
        return

    try:
        step = odb.steps[STEP_NAME]
    except KeyError:
        print("[ERROR] Step name not found in the ODB.")
        odb.close()
        return

    # Master dictionary: data[node][time][contour] = SIF_value
    data = {}

    # =====================================================================
    # 2. DATA EXTRACTION
    # =====================================================================
    for region_name, region in step.historyRegions.items():
        for ho_name, ho in region.historyOutputs.items():
            
            # NOTE: To evaluate other fracture modes, change 'K1' to 'K2' or 'K3' below.
            if 'K1' in ho_name and 'XFEM' in ho_name:
                match = re.search(r'XFEM_(\d+).*Contour\s*_?(\d+)', ho_name)
                if match:
                    node = int(match.group(1))
                    contour = int(match.group(2))
                    
                    if node not in data:
                        data[node] = {}
                        
                    for time_val, sif_val in ho.data:
                        # Round time slightly to prevent float key mismatch issues
                        t_key = round(time_val, 6)
                        if t_key not in data[node]:
                            data[node][t_key] = {}
                        data[node][t_key][contour] = sif_val

    odb.close()
    print("Data successfully extracted. Applying filtering algorithm and writing CSV...")

    # =====================================================================
    # 3. FILTERING ALGORITHM & CSV WRITING
    # =====================================================================
    with open(CSV_FILENAME, 'wb') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['XFEM_Node', 'Step_Time', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'Representative_K', 'Mesh_Status'])
        
        for node in sorted(data.keys()):
            for time_val in sorted(data[node].keys()):
                contour_values = data[node][time_val]
                
                # Prepare base row
                row = [node, time_val]
                y_list = []
                
                # Append contours C1 through C8
                for c in range(1, 9):
                    val = contour_values.get(c, None)
                    if val is not None:
                        row.append(round(val, 4))
                        y_list.append(val)
                    else:
                        row.append("N/A")
                        y_list.append(None)
                
                # Filtering Algorithm: Contiguous Moving Window
                best_mean = None
                best_error = 999.0
                best_size = 0
                
                # Iterate looking for the largest possible block (from 8 down to MIN_CONTOURS)
                for size in range(8, MIN_CONTOURS - 1, -1):
                    # Slide the window across the contours
                    for start in range(len(y_list) - size + 1):
                        window = y_list[start : start + size]
                        
                        # Skip if any data is missing in the current window
                        if None in window: 
                            continue 
                        
                        v_max = max(window)
                        v_min = min(window)
                        v_mean = sum(window) / float(size)
                        
                        if v_mean == 0: 
                            continue
                        
                        error = abs((v_max - v_min) / v_mean)
                        
                        # If it meets the tolerance, check if it is the optimal plateau
                        if error <= TOLERANCE:
                            if size > best_size:
                                best_size = size
                                best_error = error
                                best_mean = v_mean
                            elif size == best_size and error < best_error:
                                best_error = error
                                best_mean = v_mean
                    
                    # Stop searching smaller windows if a valid block was found at the current size
                    if best_size > 0:
                        break
                
                # Evaluate algorithm results and close the row
                if best_size >= MIN_CONTOURS:
                    row.append(round(best_mean, 4))
                    row.append("OK (Plateau of {} contours)".format(best_size))
                else:
                    row.append("N/A")
                    row.append("High dispersion across contours")
                    
                writer.writerow(row)

    print("Process successfully completed! File saved as: {}".format(CSV_FILENAME))

if __name__ == '__main__':
    execute_sif_post_processing()