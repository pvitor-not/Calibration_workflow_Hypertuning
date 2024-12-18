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

## FILES AND FOLDERS
# Main project directory/working Folder
# ROOT BOING
root_folder = "C:\\Users\\09008889978\\Desktop\\WorkingDir"
# ROOT LINDAURA
#root_folder = "C:\\Users\\00615045995\\PycharmProjects\\Calibration_workflow"
# ROOT JÉSSICA
#root_folder = "C:\\Users\\Jessica de Souza\\PycharmProjects\\Caliration_workflow_v2.0"
#root_folder = "C:\\Users\\45597483811\\PycharmProjects\\Calibration_workflow_v2.0_refatorado"

### STUDY CASE
case = 'reference_solution_FO_zero'

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
#Mixed_Case Lithologies
litho_class_ids = {'1':'Silt','2':'Sand','3':'Reef','4':'Carbonate'}
carbo_ids = ['3','4']
mud_text_ids = ['1']
sand_text_ids = ['2','3']
gravel_text_ids = ['4']

### END OF INPUTS ###

## Necessary folders and files definition based on main project folder 'root_folder'.
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
internal_optimization_method = 'D' #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'D'= Dijkstra
external_optimization_method = 'BF' #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'BF'= Brute Force
simulation = Objects.Simulation_parameters()
alpha = 1
beta = 1
internal_optimizer = Functions.get_optimizer(internal_optimization_method)
external_optimizer = Functions.get_optimizer(external_optimization_method)

OF = Functions.get_objective_function(OF_type)
OF.set_objective_function(case, simdir, arcfilebase, origin, well_files,extract_columns, facies_definition,simulation,
                          alpha = alpha, beta = beta,
                          litho_class_ids = litho_class_ids, carbo_ids=carbo_ids, mud_ids = mud_text_ids, sand_ids= sand_text_ids, gravel_ids=gravel_text_ids)
                                                    #PROOF PARAMETERS (OPTIONAL)
OF.set_optmization(uncertain_parameters, enginefolder, internal_optimizer, external_optimizer)

pop_size = 2 # Population size
max_iter = 2 # Maximum number of iterations (termination criterion)
best_variable, OF_value = external_optimizer.optimize(OF.compute,uncertain_parameters['min_value'].to_numpy(), uncertain_parameters['max_value'].to_numpy(),pop_size,max_iter)

Functions.plotResults(OF_type, simulation,simdir)

toc = time.time()-tac
print("Computation time calibration:", toc, "s"	)
print("The best parameter set is: ", best_variable)
print("With a OF value of: ", OF_value)