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
import Code_CW.Analysis as Analysis
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
root_folder = "C:\\Users\\10360232906\\Desktop\\Sens"
# ROOT LINDAURA
#root_folder = "C:\\Users\\00615045995\\PycharmProjects\\Calibration_workflow"
# ROOT JÉSSICA
#root_folder = "C:\\Users\\Jessica de Souza\\PycharmProjects\\Caliration_workflow_v2.0"
#root_folder = "C:\\Users\\45597483811\\PycharmProjects\\Calibration_workflow_v2.0_refatorado"

### STUDY CASE
case = 'Sand_Vol'
### DIONISOS PROJECT INPUT
# Absolute initial coordinates of the project's grid defined on Dionisos interface.
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
#M_Chaves Lithologies
litho_class_ids = {'1':'Lutites','2':'Carbo_Mud','3':'Sand','4':'Carbo_Grains','5':'Gravel','6':'Carbo_Rud'}
carbo_ids = ['2','4','6']
mud_text_ids = ['1','2']
sand_text_ids = ['3','4']
gravel_text_ids = ['5','6']
### END OF INPUTS ###

## Project Folders
## Necessary folders and files definition based on main project folder 'root_folder'.
# Project Inputs Directory
simdir = os.path.join(root_folder, case, 'input_base')
# Observed Well Files Directory
welldir = os.path.join(root_folder, case,'obs_wells')
# Uncertain Parameters File
expDesign = os.path.join(simdir,'expdesign.xlsx')
uncertain_params = os.path.join(simdir, 'uncertainties.txt')
#expDesign = ''
#uncertain_params = ''
# Facies Definition File
facies_def_file = os.path.join(simdir, 'Table_def_facies.txt')
# Well Zone Markers
well_markers_file = os.path.join(simdir, 'wellMarkers.txt')
# Reading Well Markers File
well_markers = Functions.readWellMarkers(well_markers_file)
# Reading Facies Definition Table
facies_definition = Functions.readFaciesDefinition(facies_def_file)

tac = time.time()
print("Preparation time:", tac-tic, "s")

#########################
#     FOs COMPUTING     #
#########################
# TypeOF =>'S'= Stratigraphic or 'P'= probabilistic
OF_type = 'S'
# Optimization method selection
internal_optimization_method = 'DSR' #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'D'= Dijkstra
alpha = 0
beta = 0
internal_optimizer = Functions.get_optimizer(internal_optimization_method)
# PROOF Probabilities Combining Methods (mean or belief)
evidence_method = 'mean'
total_OF_method = 'belief'
OF = Functions.get_objective_function(OF_type)
OF.internal_optimizer = internal_optimizer
simulation = Objects.Simulation_parameters()

SensAnalysis = Analysis.SensitivityAnalysis(case, OF,simdir,origin,welldir,extract_columns,facies_definition,simulation, enginefolder,
                                           litho_class_ids = litho_class_ids, carbo_ids=carbo_ids, mud_ids = mud_text_ids, sand_ids= sand_text_ids, gravel_ids=gravel_text_ids,
                                           alpha=alpha,beta=beta,internal_optimization_method=internal_optimization_method,
                                           markers_table=well_markers, expdesign =expDesign, uncertain_params = uncertain_params,
                                           evidence_combining_method=evidence_method, total_OF_computing_method=total_OF_method)
OF_Values = SensAnalysis.computeSensitivityAnalysis()
SensAnalysis.writeSensResults(OF_Values)
toc = time.time()-tac
print("Computational time spent:", toc, "s"	)
