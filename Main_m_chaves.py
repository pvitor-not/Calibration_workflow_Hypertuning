"""
Main file of the prototype of calculation of calibration workflow - 
Calibration of DionisosFlow Morro do Chaves Models
"""
##### Import libraries ####
import time
import Code_CW.Objects as Objects
import Code_CW.Functions as Functions
import Code_CW.Sim_PostProcessingAnalysis as Sim_PostProcessingAnalysis
import Code_CW.StatisticalAnalysis as Statistics
from Code_CW.Sim_PostProcessingAnalysis import RelatorioPDF
import os


#####################################
tic = time.time()

#### Inputs ####
## Dionisos Simulator (.exe) file needed to run a perform a SFM simulation ##
enginefolder = r'"C:\Program Files\Beicip\OpenFlowSuite_2023\WfEngine\plugins\fr.ifp.dionisosflow.arcades.engine.win64_15.1.0.r27540\windows\ArcaDESLauncher.exe"'

## Working path containing your calibration/analysis cases ##
root_folder = r"C:\Users\04099305201\Desktop\Modelos Dionisos"

## Calibration case folder ##

case = "FinalTest_cw"

### General Inputs ###
## Absolute initial coordinates of the project's grid defined on Dionisos interface ##
x0 = 808148
y0 = 8916612
origin = [x0,y0]

### Observed Well Inputs ###
## Name of the columns to be read on Observed Well files ##
## (Observed well must contain 4 classes of data) ##
## (Depth; Facies; Lithology; Bathymetry) ##
depth_column = 'DEPT'
facies_column = 'fac'
bathymetry_column = 'bat'
lito_column = 'lito_upscaled'
extract_columns = [depth_column, facies_column, bathymetry_column, lito_column]

## Lithology Classes Textural IDs ##
litho_class_ids = {'1':'Lutites','2':'Carbo_Mud','3':'Sand','4':'Carbo_Grains','5':'Gravel','6':'Carbo_Rud'}
carbo_ids = ['2','4','6']
mud_text_ids = ['1','2']
sand_text_ids = ['3','4']
gravel_text_ids = ['5','6']

### Project Folders ###
## Necessary folders are set based on main project folder 'root_folder' ##
simdir, welldir = Functions.set_project_folders(root_folder,case)
### Project Files ###
## Necessary files are set based on input_base path (simdir) ##
uncertain_params_file, facies_def_file, well_markers_file, zone_id_file, \
    color_ref_file, color_facies_file = Functions.set_project_files(simdir)
### Reading Necessary Files ###
## Well markers file ##
well_markers = Functions.readWellMarkers(well_markers_file)
## Facies Definition Table ##
facies_definition = Functions.readFaciesDefinition(facies_def_file)
## Observed wells files ##
well_files = Functions.getWellFiles(welldir)
## Inversible Parameters ##
uncertain_parameters = Functions.readUncertainParameters(uncertain_params_file)
## Dionisos main input file (.arc) ##
arcfilebase = Functions.getArcInput(simdir)

## End of Preparation Time ##
tac = time.time()
print("Preparation time:", tac-tic, "s")

##### Calibration Input Parameters ####
## Selecting the        ('S'= Stratigraphic Objective Function) ##
## Objective Function   ('P'= Probabilistic Objective Function) ##
OF_type = 'P'
OF = Functions.get_objective_function(OF_type)


### Objective Functions Specific Inputs ###
## SCOOF ##
# SCOOF terms weights (alpha multiplies the facies term; beta multiplies the thickness term) #
alpha = 4   # If coefs = 0, optimizer will use weighting equations;
beta = 1    # If coefs = n, for n !=0, optimizer will multiply the terms by n.

## PROOF ##
## Evidence combining methods ##
## ('mean' = average of values; 'belief' = Belief Function probability combining methods) ##
evidence_method = 'mean' # refers to the method to combine the attributes (probability evidences) included in PROOF analysis
total_OF_method = 'belief'  # refers to the method to combine well OFs (partial values) into Total OF value.


## Optimization Methods Selection ##
internal_optimization_method = 'DSR'    #'PSO'= Particle Swarm Optimization or 'GA'= Genetic Algorithm or 'D'= Dijkstra or 'DSR'= Dijkstra with with square root in the thicknesses
external_optimization_method = 'PSO'     #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'BF'= Brute Force
simulation = Objects.Simulation_parameters()
internal_optimizer = Functions.get_optimizer(internal_optimization_method)
external_optimizer = Functions.get_optimizer(external_optimization_method)

## Population Size and Number of Iterations performed by calibration method ##
pop_size = 1  # Population size
max_iter = 1  # Maximum number of Iterations (Termination criterion)

## Genetic Algorithm Aditional Parameters ##
mutation_prob = 0.05  # Mutation Probability
elit_ratio = 0.1      # Elite Ratio

# End of Inputs

### Setting Objective Function ###
OF.set_objective_function(case, simdir, arcfilebase, origin, well_files, extract_columns, facies_definition,simulation,
                          markers_table=well_markers, litho_class_ids=litho_class_ids,
                          carbo_ids=carbo_ids, mud_ids=mud_text_ids, sand_ids= sand_text_ids, gravel_ids=gravel_text_ids,
                          alpha=alpha,beta=beta,internal_optimization_method=internal_optimization_method,
                          evidence_combining_method=evidence_method, total_OF_computing_method=total_OF_method)

### Setting Optimization ###
OF.set_optmization(uncertain_parameters, enginefolder, internal_optimizer, external_optimizer)

### Computing Calibration ###
best_variable, OF_value = external_optimizer.optimize(OF.compute,uncertain_parameters['min_value'].to_numpy(), uncertain_parameters['max_value'].to_numpy(), pop_size, max_iter, args = [mutation_prob, elit_ratio])
toc = time.time()-tac
print("Computation time calibration:", toc, "s"	)
print("The best parameter set is: ", best_variable)
print("With a OF value of: ", OF_value)

### POS PROCESSING -----------------------------------------------------------------------------------------------------

# Create convergence plots on final folder
Functions.plotResults(OF_type, simulation, simdir)

#  PieCharts
domain = [x0, 8400, y0, 7500]  # [x0, xlen, y0, ylen]
print('Extrating information of wells for plotting PieCharts...')
Sim_Analysis_Prop = Sim_PostProcessingAnalysis.PropSimAnalysis(root=simdir, of_type=OF_type, well_files=well_files,
                                                               id_file=zone_id_file,
                                                               domain=domain,
                                                               extract_columns=extract_columns,
                                                               color_reference_file=color_ref_file,
                                                               facies_color=color_facies_file,
                                                               facies_def_obj=facies_definition,
                                                               well_markers_data=well_markers_file,
                                                               origin_model=origin,
                                                               arcfile=arcfilebase,
                                                               engine=enginefolder)

Sim_Analysis_Prop.wellPieCharts(Sim_Analysis_Prop.welldataprop, Sim_Analysis_Prop.id_table,
                                Sim_Analysis_Prop.sedProp)
Sim_Analysis_Prop.scatterpie(Sim_Analysis_Prop.welldataprop, Sim_Analysis_Prop.id_table,
                             Sim_Analysis_Prop.sedProp)

# Vertical well sections
Sim_Analysis_Prop.faciesExtraction(Sim_Analysis_Prop.wells_file, Sim_Analysis_Prop.cell_reference,
                                   Sim_Analysis_Prop.root_result, well_markers_file, OF_type)

# Statistical analysis
Statistics.StatisticalAnalysis(root_folder=simdir, of_type=OF_type)

# Reports
RelatorioPDF().export(root_folder=simdir, of_type=OF_type, basemodel=os.path.basename(__file__)[:-3],
                      time=round(toc/3600, ndigits=3), parameter=best_variable, value=OF_value)

