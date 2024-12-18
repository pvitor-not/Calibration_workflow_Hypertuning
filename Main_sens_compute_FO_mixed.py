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
#########################
tic = time.time()
#########################
#       INPUT DATA     #
#########################
# OpenFlow Dionisos Simulator (.exe) file needed to run a dionisos simulation.
enginefolder = '"C:\\Program Files\\Beicip\\OpenFlowSuite_2021\\WfEngine\plugins\\fr.ifp.dionisosflow.arcades.engine.win64_13.1.0.r25519\\windows\\ArcaDESLauncher.exe"'

### FILES AND FOLDERS
# Main project directory/Working Folder
# ROOT BOING
root_folder = "C:\\Users\\09008889978\\Desktop\\WorkingDir"
# ROOT LINDAURA
#root_folder = "C:\\Users\\00615045995\\PycharmProjects\\Calibration_workflow"
# ROOT JÉSSICA
#root_folder = "C:\\Users\\Jessica de Souza\\PycharmProjects\\Caliration_workflow_v2.0"
#root_folder = "C:\\Users\\45597483811\\PycharmProjects\\Calibration_workflow_v2.0_refatorado"

### STUDY CASE
case = 'Mixed_sens_compute_FO_1'

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
#M_Chaves Lithologies
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

tac = time.time()
print("Preparation time:", tac-tic, "s")

#########################
#     FOs COMPUTING     #
#########################
# TypeOF =>'S'= Stratigrafic or 'P'= probabilistic
OF_type = 'P'
# Optimization method selection
internal_optimization_method = 'DSR' #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'D'= Dijkstra
alpha = 1
beta = 1
internal_optimizer = Functions.get_optimizer(internal_optimization_method)

OF = Functions.get_objective_function(OF_type)
OF.internal_optimizer = internal_optimizer
simulation = Objects.Simulation_parameters()

SensAnalysis = Analysis.SensitivityAnalysis(case, OF,simdir,origin,welldir,extract_columns,facies_definition,simulation, enginefolder,
                                           litho_class_ids = litho_class_ids, carbo_ids = carbo_ids, clastos_grossos_ids = clastos_grossos_ids,
                                            alpha=alpha, beta =beta,
                                           reference_value=[['OutputPath', 'Output']])
OF_Values = SensAnalysis.computeSensitivityAnalysis()
SensAnalysis.writeSensResults(OF_Values)
toc = time.time()-tac
print("Computational time spent:", toc, "s"	)
