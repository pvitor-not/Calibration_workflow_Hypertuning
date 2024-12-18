# Calibration_workflow_v4.0
Repository for development of the code in python language to perform calibration  workflow - Calibration of DionisosFlow models.

This new version of the development code is under a single branch.

This repository contains the main file in python
 - Main_'case'.py - master file for each case tested/applied, as example: Main_mixed.py, Main_m_chaves.py  ...

And a folder 'Code_CW' containing python files: 
- Objects.py - with the classes of objects used for the calibration of DionisosFlow models.
- Functions.py - with the  functions/methods used for o development of calibration of DionisosFlow models.
- Dijkstra.py - with the implementation of Dijkstra's algorithm for internal optimization of SCOOF.
- ObjectiveFunction.py - with classes to set and calculation by different objective functions SCOOF, PROOF and TOTALOF.
- Optimizer.py - with classes to optimize the compute process using different methods/algorithms, as example GA, PSO, Dijkstra.
- PostProcessing.py - with classes to perform post-process calibration analysis, plot pie charts of sediment fractions, well sections, and calibration reports.
- StatisticalAnalysis.py - with classes for statistical analysis of project variables used during calibration. 
		

Cases are defined in 'calibration_workflow_cases' folder, available on 'GPER - Shared Folder'.

The folder 'cases' containing folders and files to each model case/example to be tested/calibrated with:   
- 'input_base' folder with the reference scenario simulation inputs:

        - input (folder): folder with reference scenario's input maps files.
	        
        - ouput (folder): folder with outputs of the DionisosFlow simulations.
        
        - .arc (file): reference scenario main input file with replaceable keywords (<<INV_PARAM, name = "parameter_name">>). 
 
	this folder also contains six .txt files:
	
		- 'Table_def_faceis.txt' - definition of the facies to the corresponding model

		- 'Uncertain_parameters.txt' - inversible parameters ranges definition.

		- 'wellMarkers.txt' - identifier of well limits to calibrate.
	
		- 'color_reference_sedments.txt' - RGBA reference for sedment proportion piecharts.

		- 'color_reference_facies.txt' - RGBA reference for facies to plot well section.

		- 'zone_ID.txt' - identifier of subzones for sedment proportion piecharts.
		
- 'obs_wells' folder with .las files containing well data.


**To use/test each one of the cases defined, the desired case folder must be copied to a working folder to be established on respective 'Main_File' as 'root_folder'.**

To install all the third-party Python packages, go to the terminal and run the following command:
	_pip install -r requirements_

When an additional Python package is necessary to run the code, the file requirements.txt must be updated using the 
following command:

	_pip freeze > requirements.txt_

...

To use the genetic algorithm as an optimization method, it is necessary to remove the time constraint set in the library. To do this, open the library file 'geneticalgorithm.py' and change line 66 to 'function_timeout=1000000,\'.
