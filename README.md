# Abaqus-Fatigue-Fracture-Tools
"# Abaqus Fatigue & Fracture Tools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)

A collection of advanced Python scripts developed for the Abaqus 6.14 environment (Python 2.7 compliant) to evaluate multiaxial fatigue criteria and automate the post-processing of fracture mechanics data. 

These tools were originally developed to biomechanically validate hyper-flexible 3D-printed titanium thoracic implants, but their generalized formulation makes them applicable to any structural Abaqus model (.odb) requiring advanced fatigue and fracture evaluation.

## 📦 Repository Contents

This repository contains four main automated scripts. Each script is designed to extract specific stress/strain data from the Abaqus Output Database (.odb), process it through advanced failure criteria, write the contour results back into the `.odb` for visualization, and export `.csv` reports for scientific plotting.

### 1. `Criterion_Goodman_TPMAX.py` (Maximum Principal Stress)
This script evaluates the phenomenological **Modified Goodman Criterion** by extracting the time history of the Maximum Principal Stress (sigma_1) at each integration point. 
* **Mechanics:** It calculates the eigenvalues of the symmetric stress tensor to isolate sigma_1, which is the primary driver for Mode I crack opening. 
* **Output:** Computes the stress amplitude (S_a) and mean stress (S_m), calculates the fatigue damage index, writes a new field output (`GOODMAN_INDEX`) to the ODB, and exports the top 10 most critical elements and the theoretical Goodman limit line to a CSV file.

### 2. `Criterion_Goodman_VM.py` (Signed von Mises)
A variation of the Goodman evaluation that utilizes the **Signed von Mises Equivalent Stress** instead of pure principal stresses.
* **Mechanics:** It computes the standard von Mises stress and assigns it a positive or negative sign based on the hydrostatic pressure (positive for macroscopic tension, negative for macroscopic compression). This approach is highly useful for complex multiaxial stress states where principal directions rotate during the loading cycle.
* **Output:** Identical reporting structure to the TPMAX script, providing both ODB contours and CSV data for the most critical elements.

### 3. `Criterion_Dang_Van.py` (Multiaxial Micromechanical Fatigue)
This is the most advanced fatigue script in the repository. It evaluates the micromechanical **Dang Van Multiaxial Fatigue Criterion**, which operates at the mesoscopic scale to predict infinite life limits under complex, non-proportional loading.
* **Mechanics:** 
  * Decomposes the macroscopic stress tensor into hydrostatic and deviatoric components.
  * Employs a custom, pure NumPy implementation of the **Badoiu-Clarkson geometric algorithm** to find the center of the smallest enclosing hypersphere in the 6D Mandel stress space. This identifies the microscopic residual stress tensor (rho*) that stabilizes the elastic shakedown state.
  * *Note:* This algorithm is written strictly in NumPy to maintain native compatibility with the restricted Python 2.7 environment shipped with Abaqus (which lacks SciPy).
* **Output:** Writes the `DANG_VAN_INDEX` contour to the ODB and exports the full microscopic shear stress vs. hydrostatic pressure path of the single most critical element to generate Dang Van limit diagrams.

### 4. `SIF_post_processing.py` (Stress Intensity Factor Automation)
Extracting accurate Stress Intensity Factors (SIF) from Abaqus using the Interaction Integral (Contours 1 through 8) often results in numerical dispersion. This script automates the extraction and filtering of this data for eXtended Finite Element Method (XFEM) or standard crack models.
* **Mechanics:** Extracts the SIF values (K_I, easily modifiable to K_II or K_III) across all contours. It applies a **Contiguous Moving Window Algorithm** that dynamically searches for the optimal "plateau" of stable contours.
* **Filtering:** The user defines a minimum number of grouped contours (e.g., 3) and a maximum allowed dispersion tolerance (e.g., 20%). The algorithm drops the highly distorted inner contours (close to the crack tip) and boundary-affected outer contours, returning a highly accurate, representative SIF value.
* **Output:** A structured CSV file detailing the values of C1-C8, the representative calculated K value, and the algorithmic status of the mesh plateau.

## 🚀 Usage Instructions

1. Ensure your Abaqus analysis has finished and the `.odb` file is closed in Abaqus/Viewer (to avoid read/write permission errors).
2. Open the desired Python script in any text editor.
3. Modify the **User Configuration** block at the top of the script with your specific parameters:
   * `ODB_PATH`: Name of your output database.
   * `INSTANCE_NAME`: Name of the part instance.
   * `STEP_NAME`: Name of the analysis step.
   * `CRITICAL_SET_NAME`: A predefined element set in your model (highly recommended to reduce computation time).
   * Material properties (e.g., S_ut, S_e, Dang Van constants a and b).
4. Run the script directly through the Abaqus command line:
   ```bash
   abaqus python Script_Name.py
