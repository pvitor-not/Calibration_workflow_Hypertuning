"""
Main file of the prototype of calculation of calibration workflow -
Calibration of DionisosFlow models
"""
# Import libraries
import pandas as pd
from sklearn.model_selection import ParameterGrid, ParameterSampler
from scipy.stats import uniform
import os
import time
import shutil
import Code_CW.Objects as Objects
import Code_CW.Functions as Functions

#########################
tic = time.time()
#########################
#       INPUT DATA     #
#########################
# OpenFlow Dionisos Simulator (.exe) file needed to run a dionisos simulation.
#enginefolder = '"C:\\Program Files\\Beicip\\OpenFlowSuite_2021\\WfEngine\plugins\\fr.ifp.dionisosflow.arcades.engine.win64_13.1.0.r25519\\windows\\ArcaDESLauncher.exe"'
enginefolder = r'"C:\Program Files\Beicip\OpenFlowSuite_2023\WfEngine\plugins\fr.ifp.dionisosflow.arcades.engine.win64_15.1.0.r27540\windows\ArcaDESLauncher.exe"'
## FILES AND FOLDERS
# Main project directory/working Folder
# ROOT BOING
#
#root_folder = "C:\\Users\\09008889978\\Desktop\\WorkingDir"
# ROOT LINDAURA
#root_folder = "C:\\Users\\00615045995\\PycharmProjects\\Calibration_workflow"
# ROOT JESSICA
root_folder = r"C:\Users\00584647271\Desktop\diretorio"
#root_folder = "C:\\Users\\45597483811\\PycharmProjects\\Calibration_workflow_v2.0_refatorado"

### STUDY CASE
case = "data_0"

### DIONISOS PROJECT INPUT
# Absolute initial coordinates of the project's grid defined on Dionisos interface.
# x0 = 405000
# y0 = 1055000
x0 = 808148
y0 = 8916612
origin = [x0,y0]

### GEOWELL INPUTS
# Name of the columns to be read on Observed Well files.
depth_column = 'DEPT'
facies_column = 'fac'
bathymetry_column = 'bat'
lito_column = 'lito_upscaled'
extract_columns = [depth_column, facies_column, bathymetry_column, lito_column]
del depth_column, facies_column, bathymetry_column, lito_column, x0, y0

# Observed lithology classes/ids inputs
#Mixed Case Lithologies
litho_class_ids = {'1':'Silt','2':'Sand','3':'Reef','4':'Carbonate'}
carbo_ids = ['3','4']
mud_text_ids = ['1']
sand_text_ids = ['2','3']
gravel_text_ids = ['4']

### END OF INPUTS ###
## Project Folders
## Necessary folders and files definition based on main project folder 'root_folder'.
# Project Inputs Directory
simdir = os.path.join(root_folder, case, 'input_base')
# Observed Well Files Directory
welldir = os.path.join(root_folder, case,'obs_wells')
# Uncertain Parameters File
uncertain_params_file = os.path.join(simdir, 'Uncertain_parameters.txt')
# Facies Definition File
facies_def_file = os.path.join(simdir, 'Table_def_facies.txt')
# Reading Facies Definition Table
facies_definition = Functions.readFaciesDefinition(facies_def_file)
# Getting observed wells (geowells) files from 'welldir'
well_files = Functions.getWellFiles(welldir)
# Reading Inversible Parameters
uncertain_parameters = Functions.readUncertainParameters(uncertain_params_file)
# Finding the main input file (.arc) on 'simdir' (input_base folder)
arcfilebase = Functions.getArcInput(simdir)

tac = time.time()
print ("Preparation time:", tac-tic, "s")

#########################
#   CALIBRATION LOOP    #
#########################
# COMPUTE TOTAL FO
# TypeOF =>'S'= Stratigrafic or 'P'= probabilistic
OF_type = 'S'
# Optimization method selection
internal_optimization_method = 'DSR' #'PSO'= Particle Swarm Optimization or 'GA'= Genetic Algorithm or 'SCOOFBF' = SCOOF with brute force,'D'= Dijkstra or 'DSR'= Dijkstra with with square root in the thicknesses
external_optimization_method = 'PSO' #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'BF'= Brute Force
alpha = 4
beta = 1
                                        #PROOF PARAMETER TO SET REFERENCE SOLUTION

pop_size = 5 # Population size
max_iter = 5  # Maximum number of iterations (termination criterion)

# ==========================================================
"""CALIBRATION METHOD: 'G' to Grid Search, 'R' to Random Search"""
CALIBRATION_METHOD = 'R'
# ==========================================================

if CALIBRATION_METHOD == 'G':
    print("Using the method: Grid Search")

    ## --- GRID SEARCH ---
    # param_grid = {
    #     'omega': [0, 0.25, 0.5, 0.75, 1],
    #     'phip': [0, 0.25, 0.5, 0.75, 1],
    #     'phig': [0, 0.25, 0.5, 0.75, 1],
    # }
    param_grid = {
        'omega': [0, 1],
        'phip': [0, 1],
        'phig': [0, 1],
    }
    # Generates the DataFrame with ALL grid combinations.
    df = pd.DataFrame(list(ParameterGrid(param_grid)))

elif CALIBRATION_METHOD == 'R':
    print("Using the method: Random Search")

    ## --- RANDOM SEARCH ---
    # Desired number of iterations for Random Search
    N_ITER_RANDOM = 3

    param_dist = {
        'omega': uniform(0, 1),
        'phip': uniform(0, 1),
        'phig': uniform(0, 1),
    }
    # Generates N_ITER_RANDOM random combinations.
    random_combinations = list(ParameterSampler(param_dist, n_iter=N_ITER_RANDOM, random_state=42))
    combinations = [{k: round(v, 2) for k, v in params.items()} for params in random_combinations]
    df = pd.DataFrame(combinations)

else:
    raise ValueError("CALIBRATION_METHOD it must be 'G' (Grid Search) or 'R' (Random Search).")

print(f"Total de combinações geradas: {len(df)}")

#######################################################################################################################

### SIMULATION LOOP ###
output_file = os.path.join(root_folder, case, "calibration_results.xlsx")

df['OF_value'] = None
df['time_h'] = None
for PV in uncertain_parameters["uncertain_parameter"]:
    df[PV] = None

# Standard exit path from the simulator
default_results_folder = os.path.join(simdir, "Results")

for idx, row in df.iterrows():
    omega = row['omega']
    phip = row['phip']
    phig = row['phig']

    #Destination folder for this combination
    results_folder_final = os.path.join(simdir, f"Results_{idx + 1}")

    #Cleanup: Ensures that the destination folder does not exist.
    if os.path.exists(results_folder_final):
        shutil.rmtree(results_folder_final)

    #Clears the default output folder if Dionisos does not overwrite it correctly.
    if os.path.exists(default_results_folder):
        shutil.rmtree(default_results_folder)

    simulation = Objects.Simulation_parameters()

    internal_optimizer = Functions.get_optimizer(internal_optimization_method)
    external_optimizer = Functions.get_optimizer(external_optimization_method)

    OF = Functions.get_objective_function(OF_type)
    OF.set_objective_function(case, simdir, arcfilebase, origin, well_files, extract_columns, facies_definition,
                              simulation,
                              litho_class_ids=litho_class_ids, carbo_ids=carbo_ids, mud_ids=mud_text_ids,
                              sand_ids=sand_text_ids, gravel_ids=gravel_text_ids,
                              alpha=alpha, beta=beta,
                              internal_optimization_method=internal_optimization_method)  # alpha and beta must be a fixed number or zero to use the weighting equations
    OF.set_optmization(uncertain_parameters, enginefolder, internal_optimizer, external_optimizer,
                       reference_value=[])


    tic = time.time()
    best_variable, OF_value = external_optimizer.optimize(
        OF.compute,
        uncertain_parameters['min_value'].to_numpy(),
        uncertain_parameters['max_value'].to_numpy(),
        pop_size,
        max_iter,
        [omega, phip, phig]
    )
    tac = time.time() - tic

    # ==========================================================
    # MANAGEMENT AND PROCESSING OF RESULTS
    # ==========================================================

    # The Functions.plotResults file needs to run NOW, while the folder is named 'Results'.
    # It generates the graphs within 'simdir/Results/S' or 'simdir/Results/P'.
    Functions.plotResults(OF_type, simulation, simdir)

    # EXCEL FILE CLEANUP
    if os.path.exists(default_results_folder):
        print(f"Starting Excel file cleanup for the combination {idx + 1}...")

        # Define the specific file to keep
        file_to_keep = "Uncertain_parameters_and_OF_values"

        # Cycles through the main folder ('Results') and all subfolders
        for root, dirs, files in os.walk(default_results_folder, topdown=False):
            for name in files:
                # Check if the file is an Excel file
                if name.endswith(('.xlsx', '.xls')):
                    # Delete it ONLY if it is NOT the specific file we want to keep
                    if not name.startswith(file_to_keep):
                        os.remove(os.path.join(root, name))

        print("Cleaning completed.")

    #Move and rename the 'Results' folder (now empty) to 'Results_X'.
    if os.path.exists(default_results_folder):
        shutil.move(default_results_folder, results_folder_final)
    else:
        print(
            f"Warning: Default results folder '{default_results_folder}' not found after optimization {idx + 1}. Check the simulation.")

    # ==========================================================

    #The DataFrame update continues.
    for param_name, param_value in zip(uncertain_parameters["uncertain_parameter"], best_variable):
        df.at[idx, param_name] = param_value
    df.at[idx, 'OF_value'] = OF_value
    df.at[idx, 'time_h'] = tac / 3600  # Time in hours

    print("Computation time calibration:", tac, "s")
    print("The best parameter set is: ", best_variable)
    print("With a OF value of: ", OF_value)

    #PARTIAL RESCUE
    df.insert(0, "Combination", range(1, len(df) + 1)) if "Combination" not in df.columns else None
    df.to_excel(output_file, index=False)
    print(f"\nPartial results saved in: {output_file}")

# End of loop
print(f"\nResults saved in: {output_file}")

    #######


