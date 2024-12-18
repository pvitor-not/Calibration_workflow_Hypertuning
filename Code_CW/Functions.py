"""
File of functions used for the prototype of calibration of DionisosFlow models
"""
import pandas as pd
from bs4 import BeautifulSoup as bf
import numpy as np
import os
import lasio
import itertools
import pathlib
import h5py
import shutil
import pdb
import re
from Code_CW import Dijkstra
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import coo_matrix
import Code_CW.Objects as Objects
from pyswarm import pso
import sys
import matplotlib.pyplot as plt
import Code_CW.Optimizer as Optimizer
import Code_CW.ObjectiveFunction as ObjectiveFunction
import Code_CW.xml_parser as xml_parser


def getArcInput(simdir: str):
	"""
	This function obtains the '.arc' main input file from given root simulation
	directory.
	Parameters
	----------
	simdir : str
		Folder containing temporary simulation data of the "initial" project.
	Returns
	-------
	arcfile : str
		Name of the .arc input file
	"""
	onlyfiles = [f for f in os.listdir(simdir) if os.path.isfile(os.path.join(simdir, f))]
	arcfile = [f for f in onlyfiles if pathlib.Path(f).suffix == '.arc']
	arcfile = arcfile[0]
	return arcfile

def getWellFiles(welldir: str):
	"""
	This function obtains all well files '.las' from given directory.
	Parameters
	----------
	welldir : str
		Directory containing well files.
	Returns
	-------
	onlylas : list
		List containing all well files from 'welldir'.
	"""
	onlyfiles = [f for f in os.listdir(welldir) if os.path.isfile(os.path.join(welldir, f))]
	onlylas = [os.path.join(welldir, f) for f in onlyfiles if pathlib.Path(f).suffix == '.las']
	return onlylas

def getSimData(simdir: str, file: str, origin: list, unit=1000.0):
	"""
	This function generates Dionisos_Simulation object, containing all necessary
	simulation data from the input files.
	Parameters
	----------
	simdir : str
		Folder containing temporary simulation data of the "initial" project.
	file : str
		'.arc' main input file. Can be obtained by 'get_arc_input()' function.
	origin : list
		List of two floats corresponding to the absolute X0,Y0 of the grid
	unit : float, optional.
		The default is 1000.0.
		The information obtained from mesh needs to be multiplicated for this value
		so it can match the well coordinates.
	Returns
	-------
	sim_data : object
		Object storing information on the Dionisos Simulation
	"""
	with open(os.path.join(simdir, file)) as arc:
		arc = bf(arc, 'html.parser')

	sim_data = Objects.Dionisos_Simulation(simdir, arc, origin, unit)
	return sim_data

def getMeshCoordinates(simdata:object, unit = 1000.0):
	"""
	Function designed to get mesh cells and nodes coordinates from Dionisos_Simulation object.
	Returns a list containing 'origin'/'unit' (information previously entered) and
	'xcoords'/'ycoords': lists corresponding to each mesh node (x,y) coordinate.
	This information will be used to match well positions with grid cells.
	"""
	coords = np.array([list(x) for x in simdata.mesh_data()['coords']])
	origin = simdata.origin
	xcoords = coords[:,0]
	xcoords = [x*unit for x in set(xcoords)]
	xcoords.sort()
	xcoords = xcoords[:-1]
	ycoords = coords[:,1]
	ycoords = [y*unit for y in set(ycoords)]
	ycoords.sort()
	ycoords = ycoords[:-1]
	coordinates =[origin, unit, xcoords, ycoords]
	return coordinates

def getPropFromTxtFile(wellfile, prop):
	with open(wellfile) as file:
		lines = file.read().split("\n")
		for line in lines:
			if prop in line:
				line = line.split()
				value = float(line[-1])
				break
	return value

def readWellFile(well:str, extract_columns:list):
	'''
	This function reads a .las well file based on its facies and depth column names.
	It returns a dictionary containing the following keys:
	'name': well name
	'x': x coordinate
	'y': y coordinate
	'facies': list containing facies data from well file
	'depth': list containing depth data from well files
	'''
	file = lasio.read(well)
	wellinfo = {}
	depth = extract_columns[0]
	facies = extract_columns[1]
	bathymetry = extract_columns[2]
	lithology = extract_columns[3]
	#Try to find well properties on .las file. If not found, returns a ERROR message.
	#WELL NAME
	try:
		wellinfo['name'] = file.well.WELL.value
	except AttributeError:
		print("!!ERROR: Could not find well name on .las file.!!")
	#WELL COORDINATES
	#Checking if there's _dev files on well list:
	check_file = well.replace('.las','_dev')
	if os.path.exists(check_file) == True:
		try:
			wellinfo['x'] = getPropFromTxtFile(check_file,'X-COORDINATE')
			wellinfo['y'] = getPropFromTxtFile(check_file,'Y-COORDINATE')
			print(f"Coordinates successfully obtained on _dev file for well {wellinfo['name']}")
		except AttributeError:
			print(f'!!ERROR: Could not find well {wellinfo["name"]} coordinates on _dev file.!!')
			print(f"Check coordinates data on .las/_dev file.")
			exit()
	else:
		try:
			wellinfo['x'] = file.curves.XCOORD.data
			wellinfo['y'] = file.curves.YCOORD.data
			print(f"Coordinates successfully obtained on .las file for well {wellinfo['name']}.")
		except AttributeError:
				try:
					wellinfo['x'] = file.well.XCOORD.value
					wellinfo['y'] = file.well.YCOORD.value
					print(f"Coordinates successfully obtained on .las file for well {wellinfo['name']}.")
				except AttributeError:
					print(f'!!ERROR: Theres no X,Y coordinates information on .las file.\n Check the {wellinfo["name"]} .las file.')
					exit()
	#WELL DEPTH
	try:
		wellinfo['depth'] = file.curves[f'{depth}'].data
	except KeyError:
		print(f'!!ERROR: Theres no column named {depth}.\n Check the DEPTH column name to be extracted on GeoWell {wellinfo["name"]}!!')
	#WELL FACIES
	try:
		wellinfo['facies'] = file.curves[f'{facies}'].data
	except KeyError:
		print(f'!!ERROR: Theres no column named {facies}.\n Check the FACIES column name to be extracted on GeoWell {wellinfo["name"]}!!')
	#WELL BATHYMETRY
	try:
		wellinfo['bathymetry'] = file.curves[f'{bathymetry}'].data
	except KeyError:
		print(f'!!ERROR: Theres no column named {bathymetry}.\n Check the BATHYMETRY column name to be extracted on GeoWell {wellinfo["name"]}!!')
	#WELL LITHOLOGY
	try:
		wellinfo['lithology'] = file.curves[f'{lithology}'].data
	except KeyError:
		print(f'Theres no column named {lithology}. \n Check the LITHOLOGY column name to be extracted on GeoWell {wellinfo["name"]}')
	return wellinfo

def getGeoWellData(well_files: list, meshcoordinates: list, extract_columns: list, markers_table):
	"""
	This function generates a list Geological_Well object, containing information from
	well logs.
	Parameters
	----------
	well_files : List containing well files. Can be obtained by 'get_well_files()' function.
	meshcoordinates : List containing mesh coordinates information obtained from Dionisos_Simulation object.
		Can be obtained by 'getMeshCoordinates()' function.
	extract_columns:List containing name of the columns relating facies[0] and depth[1] column on .las well file.
	Returns
	-------
	well_data : dictionary
		Dictionary relating well name (KEY) and its Geological_Well object containing
		information from its logs and relative position from simulation Grid.
	"""
	well_data = {}
	for well in well_files:
		geowell = Objects.Geological_Well(readWellFile(well,extract_columns), meshcoordinates, markers_table)
		well_data[f'{geowell.name}'] = geowell
	print('List of Geological_Wells generated containing the following wells:\n {}\n\n'.format([x for x in well_data.keys()]))
	return well_data

def getGeoWellRelativePosition(geowelldata:dict):
	'''
	This function returns the relative position of Observed Wells in relation to the Simulation Grid.
		(Needed to extract Simulated Well data correctly)
	'''
	relative_positions = {}
	for well,data in geowelldata.items():
		relative_positions[well] = data.grid_intersection()
	return relative_positions

def readFaciesDefinition(file):
	"""
	Read a string referring to the file containing the information of definition of facies
	in  a DionisosFlow simulation from their attributes
	Parameters:
	  file: facies of facies definition
	Returns a Facies_Set
	"""

	f = open(file, 'r')
	# getting the number of properties
	if f.readline()[:-1] != "Properties":
		print('Error: the file is not a file of definition of facies')
	else:
		try:
			nb_prop = int(f.readline())
		except Exception as e:
			print ('Exception:', e)

	# Getting the list of properties
	prop_list = []
	for prop in range(nb_prop):
		prop_list.append(f.readline()[:-1])

	# getting the number of facies
	if f.readline()[:-1] != "Facies":
		print('Error: the file is not a file of definition of facies')
	else:
		try:
			nb_facies = int(f.readline())
		except Exception as e:
			print ('Exception:', e)

	f.close()

	facies_table = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=(nb_prop+4))

	# Creating the dictionary of facies
	faciesDio = Objects.Facies_Set(prop_list)
	faciesDio.addfacies(facies_table)

	return faciesDio

def readUncertainParameters(file):
	"""
	Read the definition of the "inversible" parameters
	(parameters that can be modified to calibrate the model)
	Parameters:
	  file: max and min of uncertain parameters
	Returns uncertain parameters information stored in a dictionary
	"""

	f = open(file, 'r')
	# getting the number of parameters
	if f.readline()[:-1] != "#Parameters":
		print('Error: the file is not a file of uncertain parameters')
	else:
		  uncertain_parameters = pd.read_csv(file, sep ="\t+", header=None, skiprows=(1))
		  uncertain_parameters.columns = ['uncertain_parameter','min_value','max_value', 'type', 'kwargs']
	f.close()

	return uncertain_parameters

def readWellMarkers(well_markers_file):
	f = open(well_markers_file, 'r')
	# getting the number of properties
	zones = ['Well']
	if f.readline()[:-1] != "Markers":
		print('Error: the file is not a file of definition of well/zone markers ')
	else:
		try:
			nb_zones = int(f.readline())
		except Exception as e:
			print('Exception:', e)

	# Getting the list of properties
	for zone in range(nb_zones):
		zones.append(f.readline()[:-1])
	# getting the number of Wells
	if f.readline()[:-1] != "Well":
		print('Error: the file is not a file of definition of facies')
	else:
		try:
			nb_wells = int(f.readline())
		except Exception as e:
			print('Exception:', e)

	f.close()

	markers_table = pd.read_csv(well_markers_file, delim_whitespace=True, header=None, skiprows=(nb_zones + 4))
	markers_table.columns = zones

	return markers_table

def replaceInversibleParameters(replaces: list, simdir: str, arcfile: str, uncertain_parameters):
	"""
	This function reads the main input file and replaces the defined keyword
	with given values for each parameter defined in a list.
	The 'inversible' parameters must be substituted manually in the input file
	by the following keyword.
	## KEYWORD: <<INV_PARAM, name = Source1_Width>>
	The 'name' must match the parameter name on replaces list so it can be
	identified and replaced by the correct value.
	Parameters
	----------
	replaces : list
		A list of lists. Each list must have the parameter name and its corresponding
		value to be added in the input file.
	simdir : str
		Folder containing temporary simulation data of the "initial" project.
	arc_input : str
		'.arc' main input file. Can be obtained by 'get_arc_input()' function.
	Returns
	-------
	arc_input_mod: str
		The path of created '.arc' file contaning the actual parameters values
		to be used to run a Dionisos Simulation.
	"""
	arc_input = os.path.join(simdir,arcfile)
	arc_input_mod = arc_input.replace('.arc', '_mod_params.arc')
	sed_classes = sediment_classes(simdir, arc_input)
	replaces = np.array(replaces)
	indicator = 0

	with open(arc_input, 'r') as f, open(arc_input_mod, 'w') as f2:
		lines = f.readlines()
		i=0
		while i < len(lines):
			result = re.search('<<INV_PARAM(.*)>>', lines[i])
			if result != None:
				indicator = 1
				aux_txt = result.group(1).split('=')
				var = aux_txt[-1].strip()
				try:
					index = int(uncertain_parameters.index[uncertain_parameters['uncertain_parameter'] == f'{var}'].tolist()[0])
				except:
					print(f'"{var}" was defined set in .arc file as inversible but is not listed on uncertain parameters file.')
					sys.exit()
				if var not in replaces:
					print(f'Theres a inversible parameter missing on replaces list.')
					print(f"Couldn't find a input value for '{var}'.")
					f2.write(lines[i])
					i += 1
				elif uncertain_parameters.loc[index, 'type'] == 'Real':
					arcfilereplaces = uncertainScalarParser(arc_input, replaces, var, index)
					row, _ = np.where(replaces == var)
					value = arcfilereplaces[row, 1][0]
					line = re.sub('<<INV_PARAM(.*)>>', str(value), lines[i], flags=re.DOTALL)
					f2.write(line)
					i += 1
			else:
				f2.write(lines[i])
				i += 1

	for i in uncertain_parameters.index:
		type = uncertain_parameters.loc[i, 'type']
		param = uncertain_parameters.loc[i, 'uncertain_parameter']
		param_index = list(replaces[:, 0]).index(param)
		factor = replaces[param_index, 1]

		if type == 'Map':
			age = uncertain_parameters.loc[i, 'kwargs']
			uncertainMapsParser(simdir, arc_input_mod, param, age, factor)

		if type == 'Curve':
			spec = uncertain_parameters.loc[i, 'kwargs']
			uncertainCurvesParser(param, arc_input_mod, factor, spec)

		if type == 'Meta':
			age = uncertain_parameters.loc[i, 'kwargs']
			uncertainMetaParamParser(param, arc_input_mod, age, factor, sed_classes)

	if indicator == 1:
		print(f'A modified ".arc" file was created on: \n\n{simdir} \n\nwith the following values: {replaces}\n\n')
	if (indicator ==0 and any(uncertain_parameters['type']=='Map')):
		print(f'Uncertainties associated to maps were found. Modified ".h5" files were created on: \n\n {simdir}\\input \n\nwith the following values: {replaces}\n\n')
	if (indicator == 0 and any(uncertain_parameters['type'] == 'Curve')):
		print(f'Uncertainties associated to curves were found. Modified ".h5" files were created on: \n\n {simdir}\\input \n\nwith the following values: {replaces}\n\n')
	if (indicator == 0 and any(uncertain_parameters['type'] == 'Meta')):
		print(f'Uncertainties associated to metaparameters were found. Modified ".h5" files were created on: \n\n {simdir}\\input \n\nwith the following values: {replaces}\n\n')
	if (indicator == 0 and (all(uncertain_parameters['type']!= 'Map') and all(uncertain_parameters['type']!='Curve') and all(uncertain_parameters['type']!='Meta'))):
		print('No inversible parameters found in the .arc input file represented by the correct keyword.')
		print('\nInsert inversible parameters using: \nKEYWORD: <<INV_PARAM, name = "Parameter ID">>')
		print('\n"Parameter ID" on ".arc" file must match "replaces" list given ID')
		print('Modified ".arc" input file generated is equal to ".arc" base file.')
		exit()
	for i in range(len(replaces)-1):
		age = uncertain_parameters.loc[i,"kwargs"]
		if age != 0:
			replaces[i,0] = f'{uncertain_parameters.loc[i,"uncertain_parameter"]}{age}'
	arc_mod = arcfile.replace('.arc', '_mod_params.arc')
	return arc_mod,replaces

def getMetaParamData(arc_dict, param, age):
	if ('Ratio' in param or 'ratio' in param or 'proportion' in param or 'Proportion' in param or 'prop' in param or 'Prop' in param):
		sed_source_data = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s']['stratigraphy-boundary-data']['boundary-supplies']['age'],'')
		nb_sources = int(sed_source_data[0]['nb-sources'])
		for k in range(len(sed_source_data)):
			time = abs(float(sed_source_data[k]['age']))
			if time == abs(float(age)):
				age_idx = k
				for i in range(0, nb_sources):
					source_id = str(i + 1)
					if f'S{source_id}' in param or f'Source{source_id}' in param or f'source{source_id}' in param:
						source_idx = i
						df_data = pd.DataFrame(data = sed_source_data[k]['source'][i]['concentration-of-lithology'], columns = ['ratio'])
						df_data= pd.DataFrame(df_data.loc[:,'ratio'].astype(float))
						to_drop = df_data.index[df_data.loc[:,'ratio'] == 0].tolist()
						df_data_master_slave = df_data.drop(df_data.index[to_drop])
		kwargs = [age_idx, source_idx]
		return df_data, df_data_master_slave, kwargs

def computeMetaParamValues(df_data, master_slave_data, sed_classes, param, pct):
	master = [key for key,item in sed_classes.items() if item in param][0]
	master_slave_data.loc[:, 'relation'] = np.zeros(len(master_slave_data.ratio))
	if len(master_slave_data) == 2:
		for i in master_slave_data.index:
			if i == master:
				df_data.loc[i,'ratio'] = pct
			else:
				df_data.iloc[i, 0] = 1 - pct

	if len(master_slave_data) >2:
		res = 1 - master_slave_data.loc[master, 'ratio']
		new_res = 1 - pct
		for i in master_slave_data.index:
			if i != master:
				master_slave_data.loc[i,'relation'] = master_slave_data.loc[i,'ratio']/res
				df_data.loc[i,'ratio'] = new_res * master_slave_data.loc[i, 'relation']
			else:
				master_slave_data.loc[i, 'relation'] = 1
				df_data.loc[i,'ratio'] = pct

	return df_data

def replaceMetaParam(arc_input, new_meta_data, param, kwargs = []):
	new_xml_file = xml_parser.xml_to_dict(arc_input)
	new_meta_data =[f'{x:.8f}' for x in new_meta_data.ratio]
	if ('Ratio' in param or 'ratio' in param or 'proportion' in param or 'Proportion' in param):
		age_idx, source_idx = kwargs
		new_xml_file['case']['arca-d-e-s']['stratigraphy-boundary-data']['boundary-supplies']['age'][int(f'{age_idx}')]['source'][int(f'{source_idx}')]['concentration-of-lithology'] = new_meta_data
	new_xml_file = xml_parser.dict_to_xml_in_string(new_xml_file)
	with open(arc_input,'w') as xml:
		xml.write(new_xml_file)
		xml.close()

def uncertainMetaParamParser(param, arc_input, age, pct, sed_classes):
	arc_dict = xml_parser.xml_to_dict(arc_input)
	pct = float(pct)/100
	df_data, master_slave_data, kwargs = getMetaParamData(arc_dict, param, age)
	new_meta_data = computeMetaParamValues(df_data, master_slave_data, sed_classes, param, pct)
	replaceMetaParam(arc_input, new_meta_data,param, kwargs)
	return df_data

def uncertainCurvesParser(param, arc_input_mod, factor, spec):
	xml_curve,param_id = getParameterCurve(param, arc_input_mod, spec)
	writeInversibleCurve(param, param_id, arc_input_mod, xml_curve, factor, spec)

def getParameterCurve(param, arc_input, spec):
	arc_dict = xml_parser.xml_to_dict(arc_input)
	param_id = ''
	if ('Eustasy' or 'eustasy') in param:
		xml_curve = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s']['eustasy-data'], ['eustasy-time-data'])
	elif ('Prod' or 'prod') in param:
		litho_data = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s'],['lithology'])

		for litho in litho_data:
			if (param in litho['litho-name'] or litho['litho-name'] in param):
				param_id = litho['identifier']
				continue

		if spec =='vsBat':
			litho_prod_data = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s']['litho-production']['production-time-bathy-dependency'],['prod-data-for-litho'])
			xml_curve = litho_prod_data[int(param_id)]['bathy-production-for-litho']['prod-bathy-data']

		elif spec =='vsTime':
			litho_prod_data = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s']['litho-production']['production-time-bathy-dependency'],['prod-data-for-litho'])
			xml_curve = litho_prod_data[int(param_id)]['time-production-for-litho']['prod-time-data']

	elif (('Supply' or 'supply') in param or ('Volume' or 'volume') in param):
		sed_source_data = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s']['stratigraphy-boundary-data']['boundary-supplies']['age'], '')
		nb_sources = int(sed_source_data[0]['nb-sources'])
		xml_curve = []
		for k in range(len(sed_source_data)):
			for i in range(0, nb_sources):
				source_id = str(i + 1)
				if f'S{source_id}' in param or f'Source{source_id}' in param or f'source{source_id}' in param:
					xml_curve.append(float(sed_source_data[k]['source'][i]['height']))
					param_id = i

	return xml_curve, param_id

def writeInversibleCurve(param, param_id, arc_input, xml_curve, factor, spec):
	new_xml_file = xml_parser.xml_to_dict(arc_input)
	new_xml_curve = []
	if ('Eustasy' or 'eustasy') in param:
		for item in xml_curve:
			item_dict = {}
			item_dict['eustasy'] = str(float(item['eustasy'])*float(factor))
			item_dict['time'] = item['time']
			new_xml_curve.append(item_dict)
		new_xml_file['case']['arca-d-e-s']['eustasy-data']['eustasy-time-data'] = new_xml_curve

	if ('Prod' or 'prod') in param:
		if spec =='vsBat':
			for item in xml_curve:
				item_dict = {}
				item_dict['bathy'] = item['bathy']
				item_dict['prod'] = str(float(item['prod'])*float(factor))
				new_xml_curve.append(item_dict)
				new_xml_file['case']['arca-d-e-s']['litho-production']['production-time-bathy-dependency']['prod-data-for-litho'][int(param_id)]['bathy-production-for-litho']['prod-bathy-data'] = new_xml_curve
		elif spec =='vsTime':
			for item in xml_curve:
				item_dict = {}
				item_dict['time'] = item['time']
				item_dict['prod'] = str(float(item['prod'])*float(factor))
				new_xml_curve.append(item_dict)
				new_xml_file['case']['arca-d-e-s']['litho-production']['production-time-bathy-dependency']['prod-data-for-litho'][int(param_id)]['time-production-for-litho']['prod-time-data'] = new_xml_curve
	elif (('Supply' or 'supply') in param or ('Volume' or 'volume') in param):
		for item in range(0,len(xml_curve)):
			new_xml_file['case']['arca-d-e-s']['stratigraphy-boundary-data']['boundary-supplies']['age'][item]['source'][param_id]['height'] = str(float(xml_curve[item])*float(factor))
	new_xml_file = xml_parser.dict_to_xml_in_string(new_xml_file)
	with open(arc_input,'w') as xml:
		xml.write(new_xml_file)
		xml.close()

def uncertainMapsParser(simdir, arcfile,param, age, factor):
	param, map, h5file, h5adress, xml_file = getParameterMap(param, age, simdir, arcfile)
	h5_newfile = writeInversibleMap(simdir, map, h5file, h5adress, factor)
	replaceUncertainMapPath(simdir, age, xml_file, h5adress, h5_newfile)

def writeInversibleMap(simdir, map, h5file, h5adress, factor):
	new_map = map[:]*float(factor)
	h5filenew = h5file.replace('.h5','_mod.h5')
	shutil.copyfile(f'{simdir}\\{h5file}', f'{simdir}\\{h5filenew}')
	with h5py.File(f'{simdir}\\{h5filenew}', 'r+') as original_map:
		original_map[f'{h5adress}'][...] = new_map
		original_map.close()
	return h5filenew

def replaceUncertainMapPath(simdir, param_age, xml_file, h5adress, newh5path):
	new_xml_file = xml_parser.xml_to_dict(f'{simdir}\\{xml_file}')
	for i in range (0,len(new_xml_file['Xdmf']['Domain']['Grid']['Grid'])):
		time = new_xml_file['Xdmf']['Domain']['Grid']['Grid'][i]['Time']['@Value']
		if abs(float(time)) == abs(float(param_age)):
			new_xml_file['Xdmf']['Domain']['Grid']['Grid'][i]['Attribute']['DataItem']['#text'] = f'{newh5path}:{h5adress}'
			continue
	new_xml_file = xml_parser.dict_to_xml_in_string(new_xml_file)
	newxml_filename = xml_file.replace('.xmf','_mod.xmf')
	with open(f'{simdir}\\{newxml_filename}','w') as xml:
		xml.write(new_xml_file)
		xml.close()
	return newxml_filename

def getParameterMap(param, param_age, simdir, arc_input):
	arc_dict = xml_parser.xml_to_dict(arc_input)
	if ('Subs' or 'subs') in param:
		xml_map = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s']['substratum-data']['subsidence-data'], ['filename'])
		subsidence_xml_adress = arc_dict['case']['arca-d-e-s']['substratum-data']['subsidence-data']['filename']
		arc_dict['case']['arca-d-e-s']['substratum-data']['subsidence-data']['filename'] = subsidence_xml_adress.replace('.xmf','_mod.xmf')
		new_arc_xml = xml_parser.dict_to_xml_in_string(arc_dict)
		with open(arc_input, 'w') as arc:
			arc.write(new_arc_xml)
			arc.close()


	map_xml = xml_parser.xml_to_dict(f'{simdir}\\{xml_map}')
	h5file = 'input\\' + xml_parser.deep_keys(map_xml['Xdmf']['Domain']['Geometry']['DataItem'][0], ['#text']).split(':')[0]

	for grid in xml_parser.deep_keys(map_xml['Xdmf']['Domain']['Grid'],['Grid']):
		time = grid['Time']['@Value']
		if abs(float(time)) == abs(float(param_age)):
			adress = grid['Attribute']['DataItem']['#text'].split(':')[-1]
			continue

	map = xml_parser.reader_h5(f'{simdir}\\{h5file}', adress)

	return param, map, h5file, adress, xml_map

def uncertainScalarParser(arc_input, replaces,param,index):
	arc_dict = xml_parser.xml_to_dict(arc_input)
	newreplaces = replaces.copy()
	if (('Supply' or 'supply') in param or ('Volume' or 'volume') in param):
		volume = float(replaces[index,1])
		sed_source_data = xml_parser.deep_keys(arc_dict['case']['arca-d-e-s']['stratigraphy-boundary-data']['boundary-supplies']['age'],'')
		nb_sources = int(sed_source_data[0]['nb-sources'])
		sed_source_widths = {}
		for i in range (0, nb_sources):
			source_id = str(i+1)
			sed_source_widths[source_id] = float(sed_source_data[0]['source'][i]['width'])
			if f'S{source_id}' in param or f'Source{source_id}' in param or f'source{source_id}' in param:
				height = (volume*2000)/ sed_source_widths[source_id]
				newreplaces[index,1] = str(height)
	return newreplaces

def toWellExtraction(simdir, file, origin, extract_columns):
	'''
	This function access Dionisos_Simuation object and extracts the attribute 'sim_data_to_extract()'.
	Returns a dictionary containing input/output files information needed to extract simulated well results.
	'''
	simdata_i = getSimData(simdir, file, origin)
	sim_info_to_extract = simdata_i.sim_data_to_extract()
	sim_info_to_extract['extract_columns'] = extract_columns
	return sim_info_to_extract

def getSimWellData(well_files:str, simdir:str, relative_positions: dict, sim_info_to_extract: dict, faciesdef):
	"""
	This function generates Simulated_Well object, containing information
	extracted from the results of a Dionisos simulation.
	Parameters
	----------
	well_files : List
		List containing well files. Can be obtained by 'get_well_files()' function.
	simdata : object
		Object containing information of a Dionisos simulation.
	geowelldata : dict
		Dictionary containing all Geological_Well objects.
	Returns
	-------
	well_data : dictionary
		Dictionary relating well name (KEY) and its Simulated_Well object containing
		information the results from Dionisos Simulation for cells intercepted by
		that well.
	"""
	well_data = {}
	for well in well_files:
		simwell = Objects.Simulated_Well(readWellFile(well, sim_info_to_extract['extract_columns']),relative_positions, sim_info_to_extract, faciesdef)
		well_data[f'{simwell.name}'] = simwell
	print('List of Simulated_Wells generated containing the following wells:\n {}\n\n'.format([x for x in well_data.keys()]))
	return well_data

def faciesDistances(values, faciesdef):
	"""
	Function to compute the distance of a simulated cell to each geological facies
	and determine the nature of the cell facies based on cell properties

	Parameters:
		values: Series of properties
		faciesdef: Facies_Set
	Returns a Series of distances to each facies
		and the name of the facies corresponding to the property values (string)
		returns 'Undetermined' if no facies corresponds to the cell properties
	"""

	list_dist = []
	cellf = 'Undetermined'
	for j, facies in faciesdef.faciesTab.iterrows():
		test = True
		dist = 0
		for prop in values.index:
			cmin = prop+'_min'
			cmax = prop+'_max'
			var_dist = 0.
			if cmin not in faciesdef.faciesTab.columns:
				continue # property not used to define facies
			elif facies[cmin] != None and facies[cmin] > values[prop]:
				test = False # Not the right facies => compute distance
				if facies[cmin] != None:
					var_dist = (facies[cmin]-values[prop])**2
			elif values[prop]>facies[cmax] and facies[cmax] != None:
				test = False
				if facies[cmax] != None:
					var_dist = (facies[cmax]-values[prop])**2
			dist += var_dist

		if test:
			cellf = facies.name
			list_dist.append(0)
		if not test:
			list_dist.append(dist**(1/2.))
	return cellf, list_dist

def runDionisosSim(root: str, engine: str, inputfile: str):
	"""
	Function to run Dionisos Sim Engine without the OpenFlow Interface.
	Output Results will be stored in ProjectRootFolder path.
	----------
	Parameters
	engine : Dionisos Simulation Engine (.exe) file path
	root : Folder containing temporary data of the project
	inputfile : .arc project input file
	"""
	root = root.replace('C:\\','"C:\\')
	inputfile = inputfile.replace('.arc','.arc"')
	inputfile = f'\\{inputfile}'
	CommandOptions = ' -t ="quiet" -o ' + root + 'arcades-output.log" --mpirun "C:\\Program Files\\Beicip\\OpenFlowSuite_2021\\WfEngine\\plugins\\fr.ifp.openflowsuite.intelmpi.win64_2019.3.203\\IntelMPI\\intel64\\bin\\mpiexec.exe" --mpirun-options=--delegate'
	ShellCommand = 'cmd /c "cd ' + root + '" && ' + engine + ' -n 1 '  + root + inputfile + CommandOptions + '"'
	os.system(ShellCommand)

def depth_to_thickness(sequence):
	"""
	Function to transform depth (MD) in .las files into thickness

	Parameters:
		depth: List of depth (MD)
	Returns the transformed thickness for the data without the final duplication
	"""

	thickness = []
	if (np.sort(sequence['depth']) == sequence['depth']).all() or (-np.sort(-sequence['depth']) == sequence['depth']).all():
		for layer in range(len(sequence['depth']) - 1):
			thickness.append(abs(
				round(abs(sequence['depth'].iloc[layer] - sequence['depth'].iloc[layer + 1]), 10)))
		sequence.drop(sequence.tail(1).index, inplace=True)
		sequence['Thickness'] = thickness
	else:
		sequence['Thickness'] = sequence['depth']
		print("Warning: data is already setted to Thickness.")
	return (sequence)

def sediment_classes(simdir,file):
	"""
	Returns a list containing sediment classes present on simulation.
	The sediment position in list corresponds to its identifier in .arc
	input file.
	"""
	with open(os.path.join(simdir, file)) as arc:
		arc = bf(arc, 'html.parser')
	litho = list(arc.find_all('lithology'))
	nfacies = len(litho)
	litholist = {}
	for x in litho:
		identifier = x.identifier
		identifier = int(identifier.text)
		sclass = x.find('litho-name')
		sclass = sclass.text
		litholist[identifier] = sclass
	return litholist

def bruteForce(func, lb, ub, samples_number, args):
	"""
	 Method to compute the total OF with a homogeneous uncertain parameters sampling
	Parameters:
	  func: arguments passed to the objective function
	  lb: lower bounds of the design variable
	  ub: upper bounds of the design variable
	  samples_number: number of samples to generate
	  args: additional arguments passed to OF
	Returns the minimum total OF value and the corresponding parameters found
	"""

	if lb.size != ub.size:
		sys.exit("Error: Lower and upper bounds sizes must be the same!")

	if samples_number < 2:
		sys.exit("Error: Number of samples must be equal or greater than 2!")

	uncertain_parameters_sample = []
	for parameter in range(lb.size):
		uncertain_parameters_sample.append(np.linspace(lb[parameter],ub[parameter],samples_number))

	uncertain_parameters_combinations = itertools.product(*tuple(uncertain_parameters_sample))

	total_OF_Values = []
	for combination in uncertain_parameters_combinations:
		total_OF_Values.append(func(np.asarray(combination),*args))

	uncertain_parameters_list = list(itertools.product(*tuple(uncertain_parameters_sample)))

	return (uncertain_parameters_list[total_OF_Values.index(min(total_OF_Values))],min(total_OF_Values))

def get_objective_function(OF_type):
	"""
	Method to create the ObjectiveFunction object based on the OF type chosen by the user

	Parameters:
		OF_type: OF type chosen by the user
	Return the object created
	"""
	if OF_type == "S":
		return ObjectiveFunction.Stratigraphic_OF()
	elif OF_type == "P":
		return ObjectiveFunction.Probabilistic_OF()
	else:
		sys.exit('Error: OF type must be S or P')

def get_optimizer(optimizer_type):
	"""
	Method to create the Optimizer object based on the optimizer type chosen by the user

	Parameters:
		optimizer_type: Optimizer type chosen by the user
	Return the object created
	"""
	if optimizer_type == "AG":
		return Optimizer.Genetic_Algorithm()
	elif optimizer_type == "PSO":
		return Optimizer.Particle_Swarm_Optimization()
	elif optimizer_type == "SCOOFBF":
		return Optimizer.SCOOF_Brute_Force()
	elif optimizer_type == "D":
		return Optimizer.Dijkstra()
	elif optimizer_type == "DSR":
		return Optimizer.Dijkstra_With_Square_Root_Thickness()
	elif optimizer_type == "BF":
		return Optimizer.Brute_Force()
	else:
		sys.exit('Error: OF type must be AG, PSO, D, DSR or BF')

def plotResults(OF_type, simulation,simdir):
	"""
	Method to plot the results

	Parameters:
	  simulation: object with all output data
	Returns the data plotted
	"""
	results_output = f'{simdir}\\Results\\{OF_type}'
	if os.path.isdir(results_output) == False:
		os.makedirs(results_output)

	np = len(simulation.alloutputdata[0])
	parameters = simulation.header[:np-1]
	lastsim = simulation.alloutputdata[len(simulation.alloutputdata)-1][len(simulation.header)-1]

	for count, parameter in enumerate(parameters):
		plt.clf()
		if parameter == 'OF Value':
			df_plot = pd.DataFrame({'x': [int(item[np - 1]) for item in simulation.alloutputdata],
									'y':[float(item[count]) for item in simulation.alloutputdata],
									'c': [int(item[np - 1]) for item in simulation.alloutputdata]})
			plt.scatter(df_plot.x, df_plot.y, s=50, c=df_plot.c, cmap = 'plasma', edgecolors = 'k')
			plt.colorbar()
			plt.xlabel('Simulations')
			plt.xlim([0,lastsim])
			plt.ylabel('OF Value')
			plt.savefig(f'{results_output}\\_@{parameter}.png')
		elif parameter =='OutputPath':
			pass
		else:
			df_plot = pd.DataFrame({'x': [float(item[count]) for item in simulation.alloutputdata],
									'y': [float(item[np - 2]) for item in simulation.alloutputdata],
									'c': [x[np - 1] for x in simulation.alloutputdata]})
			plt.scatter(df_plot.x, df_plot.y, s=50, c=df_plot.c, cmap = 'plasma', edgecolors = 'k')
			plt.colorbar()
			plt.xlabel(str(parameter))
			plt.ylabel('OF Value')
			plt.savefig(f'{results_output}\\_{parameter}.png')
		plt.clf()

	return

def set_project_folders(root_folder, case):
	simdir = os.path.join(root_folder, case, 'input_base')
	welldir = os.path.join(root_folder, case, 'obs_wells')
	return simdir, welldir

def set_project_files(simdir):
	uncertain_params_file = os.path.join(simdir, 'Uncertain_parameters.txt')
	facies_def_file = os.path.join(simdir, 'Table_def_facies.txt')
	well_markers_file = os.path.join(simdir, 'wellMarkers.txt')
	zone_id_file = os.path.join(simdir, 'zone_ID.txt')
	color_ref_file = os.path.join(simdir, 'color_reference.txt')
	color_ref_file_facies = os.path.join(simdir, 'color_reference_facies.txt')
	return uncertain_params_file, facies_def_file, well_markers_file, zone_id_file, color_ref_file, color_ref_file_facies
