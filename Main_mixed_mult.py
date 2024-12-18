"""
Main file of the prototype of calculation of calibration workflow - 
Calibration of DionisosFlow models
"""
# Import libraries
import pandas as pd
import os
import time
import Code_CW.Objects as Objects
import Code_CW.Functions as Functions
from prettytable import PrettyTable
import statistics
import matplotlib.pyplot as plt
#########################
tic = time.time()
#########################
#       INPUT DATA     #
#########################
# OpenFlow Dionisos Simulator (.exe) file needed to run a dionisos simulation.
enginefolder = '"C:\\Program Files\\Beicip\\OpenFlowSuite_2021\\WfEngine\plugins\\fr.ifp.dionisosflow.arcades.engine.win64_13.1.0.r25519\\windows\\ArcaDESLauncher.exe"'

## FILES AND FOLDERS
# Main project directory/working Folder
# ROOT BOING
#
#root_folder = "C:\\Users\\09008889978\\Desktop\\WorkingDir"
# ROOT LINDAURA
#root_folder = "C:\\Users\\00615045995\\PycharmProjects\\Calibration_workflow"
# ROOT JESSICA
root_folder = "C:\\Users\\45597483811\\PycharmProjects\\Calibration_workflow_v3.0\\cases"
#root_folder = "C:\\Users\\45597483811\\PycharmProjects\\Calibration_workflow_v2.0_refatorado"

### STUDY CASE
case = 'reference_solution'

### DIONISOS PROJECT INPUT
# Absolute initial coordinates of the project's grid defined on Dionisos interface.
x0 = 405000
y0 = 1055000
origin = [x0,y0]

### GEOWELL INPUTS
# Name of the columns to be read on Observed Well files.
depth_column = 'MD'
facies_column = 'Facies_Petro'
bathymetry_column = 'Bat_Code'
lito_column = 'Litho_Code'
extract_columns = [depth_column, facies_column, bathymetry_column, lito_column]
del depth_column, facies_column, bathymetry_column, lito_column, x0, y0

# Observed lithology classes/ids inputs
#Mixed Case Lithologies
litho_class_ids = {'1':'Silt','2':'Sand','3':'Reef','4':'Carbonate'}
carbo_ids = ['3','4']
clastos_grossos_ids = ['2','3','4']

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
facies_def_file = os.path.join(simdir, 'Table_def_faceis.txt')
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
internal_optimization_method = 'D' #'PSO'= Particle Swarm Optimization or 'GA'= Genetic Algorithm or 'SCOOFBF' = SCOOF with brute force,'D'= Dijkstra or 'DSR'= Dijkstra with with square root in the thicknesses
external_optimization_method = 'PSO' #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'BF'= Brute Force
simulation = Objects.Simulation_parameters()

internal_optimizer = Functions.get_optimizer(internal_optimization_method)
external_optimizer = Functions.get_optimizer(external_optimization_method)

OF = Functions.get_objective_function(OF_type)
OF.set_objective_function(case, simdir, arcfilebase, origin, well_files,extract_columns, facies_definition,simulation,
                          litho_class_ids=litho_class_ids, carbo_ids=carbo_ids, clastos_grossos_ids=clastos_grossos_ids,alpha=0,beta=0,internal_optimization_method=internal_optimization_method) #alpha and beta must be a fixed number or zero to use the weighting equations
                                                  #PROOF PARAMETERS (OPTIONAL)
#REFERENCE VALUES (must be given as a list of lists contaning: ['parameter name', parametervalue]):
#1) [['Sand_KwMarine', 130.5]] | mixed_FO_Coef_Dif Case
#2) [['S1_Sed_Volume', 2240]]   | mixed_FO_Sed_Vol Case
#3) [['S1_70_Height', 12], ['S2_70_Height', 28], ['S2_65_Height',10], ['S1_28_Height',20]]   | reference_solution Case
OF.set_optmization(uncertain_parameters, enginefolder, internal_optimizer, external_optimizer,
                  reference_value =  [['S1_70_Height', 12], ['S2_70_Height', 28], ['S2_65_Height',10], ['S1_28_Height',20]])
                                        #PROOF PARAMETER TO SET REFERENCE SOLUTION

pop_size = 2 # Population size
max_iter = 2 # Maximum number of iterations (termination criterion)

total_runs = 2
resultHeader = ["Run", "Computation time calibration (s)", "Best parameter set", "OF value"]

resultFileData = PrettyTable(resultHeader)
resultFileData.padding_width = 1

allOFValues = []

resultdir = os.path.join(root_folder, "results")
if not os.path.exists(resultdir):
    os.mkdir(resultdir)

for run_number in range(total_runs):
    simResultdir = os.path.join(resultdir, str(run_number))
    if not os.path.exists(simResultdir):
        os.mkdir(simResultdir)
    initialTime = time.time()
    best_variable, OF_value = external_optimizer.optimize(OF.compute,uncertain_parameters['min_value'].to_numpy(), uncertain_parameters['max_value'].to_numpy(), pop_size, max_iter)
    allOFValues.append(OF_value)

    Functions.plotResults(OF_type, simulation,simResultdir)

    toc = time.time() - initialTime

    resultFileData.add_row([run_number, toc, best_variable, OF_value])

    with open(f'{resultdir}\\All_runs_results.txt', 'w') as resultFile:
        resultFile.write(str(resultFileData))

    print("Computation time calibration:", toc, "s")
    print("The best parameter set is: ", best_variable)
    print("With a OF value of: ", OF_value)

with open(f'{resultdir}\\Statistic_results.txt', 'w') as f:
    print("OF values arithmetic mean: ", statistics.mean(allOFValues), file=f)
    print("OF values sample variance:", statistics.variance(allOFValues), file=f)
    print("OF values population variance: ", statistics.pvariance(allOFValues), file=f)

plt.title('OF values histogram')
plt.xlabel('OF values')
plt.ylabel('Absolute frequency')
plt.hist(allOFValues, 10, rwidth=0.9)
plt.savefig(f'{resultdir}\\histogram.png')
plt.show()