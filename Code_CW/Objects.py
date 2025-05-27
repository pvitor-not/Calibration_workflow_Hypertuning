"""
Classes of objects used for the calibration of DionisosFlow models
	Dionisos_Simulation: object storing information on the Dionisos Simulation
	Calibration_Set: set of simulated pseudo-wells and corresponding geological wells , at the desired resolution, used as calibration data
	Geological_Well: geological information (thickness, facies...) coming from well logs
	Simulated_Well: data from a pseudo-well extracted from a DionisosFlow simulation at a well location
"""

#Import libraries
import pandas as pd
from pandas import DataFrame as df
from bs4 import BeautifulSoup as bf
import os
import lasio
import h5py
import numpy as np
import Code_CW.Functions as Functions
from prettytable import PrettyTable
import matplotlib.pyplot as plt
import matplotlib

class Facies_Set:
	def __init__(self, list_of_prop):
		"""
		Initialize the facies_set object.
		Parameters:
		  list_of_prop : list of prop used to define facies
		"""
		self.nb_prop = len(list_of_prop)
		self.prop_list = list_of_prop
		self.faciesTab = None # dataframe of normalized facies characteristics
		self.norms = None # dictionary of normalization parameters
		
	def addfacies(self, table):
		"""
		Add the facies to the dict.
		Parameters:
		  table : dataframe of facies definition
		Returns a dataframe of Facies (normalized criteria) (faciesTab)
		Returns a dictionary of normlization parameters (norms)
		"""
		norms = {}
		new_table = pd.DataFrame()
		for prop in range(self.nb_prop):
			# print table
			pos = 1+prop*2
			mini = np.nanmin(table[table.columns[[pos,pos+1]]])
			maxi = np.nanmax(table[table.columns[[pos,pos+1]]])
			factor = maxi-mini
			new_table[self.prop_list[prop] + '_min'] = (-mini + table[pos]) / factor
			new_table[self.prop_list[prop] + '_max'] = (-mini + table[pos+1]) / factor
			norms.update({self.prop_list[prop]: (mini, factor)})
		
		new_table = new_table.set_index(table[0])
		self.norms = norms
		self.faciesTab = new_table


class Dionisos_Simulation:
	"""
	Summarize the information of a DionisosFlow simulation
	"""
	def __init__(self, simdir: str, simfile: str, origin: list, unit: float):
		"""
        Parameters
        ----------
        simdir : str
            Folder containing simulation data of the "initial" project and other project inputs.
        simfile : str
            '.arc' main input file. Can be obtained by 'get_arc_input()' function.
        origin : list
            List of two floats corresponding to the absolute X0,Y0 of the grid
        unit : float
            The default is 1000.0.
            The information obtained from mesh needs to be multiplicated for this value
            so it can match the well coordinates.
        """
		self.file = simfile
		self.simdir = simdir
		self.origin = origin
		self.x0 = origin[0]
		self.y0 = origin[1]
		self.unit = unit
		self.dx = None # TO DO
		self.dy = None # TO DO
		self.Lx = None # TO DO
		self.Ly = None # TO DO

	def final_age(self):
		"""
        Gets simulation final timestep age.
        """
		endtime = list(self.file.find_all('end-time'))
		endtime = float(endtime[0].text)
		return endtime

	def initial_age(self):
		"""
        Gets simulation intial timestep age.
        """
		init = list(self.file.find_all('initial-age'))
		init = float(init[0].text)
		return init

	def timesteps(self):
		"""
        Returns a list containing all simulation timesteps.
        """
		timesteps = []
		timesteps.append(float(self.initial_age()))
		intermediate = list(self.file.find_all('intermediate-age'))
		timesteps += [float(x.text) for x in intermediate]
		timesteps.append(float(self.final_age()))
		return timesteps

	def mesh_path(self):
		"""
        Returns the path of '.ixm' file describing how mesh information is
        stored in 'mesh'.h5 file.
        """
		path = self.file.mesh.text.strip()
		return path

	def input_maps(self):
		"""
        Returns a list containing all map input data files.
        """
		inputs = [x.text for x in self.file.find_all('filename')]
		return inputs

	def output_folder(self):
		"""
        Returns the name of the folder containing output data.
        """
		outputpath = list(self.file.find_all('output-path'))
		outputpath = outputpath[0].text.strip()
		return outputpath

	def base_file_name(self):
		"""
        Returns the default name for output files.
        """
		base = list(self.file.find_all('base-file-name'))
		base = base[0].text.strip()
		return base

	def output_variables(self):
		"""
        Return a list contanining all output variables.
        """
		outputvar = [var.text for var in self.file.find_all('variable')]
		return outputvar

	def sediment_classes(self):
		"""
        Returns a list containing sediment classes present on simulation.
        The sediment position in list corresponds to its identifier in .arc
        input file.
        """
		litho = list(self.file.find_all('lithology'))
		nfacies = len(litho)
		litholist = list(np.zeros(nfacies))
		for x in litho:
			identifier = x.identifier
			identifier = int(identifier.text)
			sclass = x.find('litho-name')
			sclass = sclass.text
			litholist[identifier] = sclass
		return litholist

	def n_sed_classes(self):
		"""
        Returns the number of sediment classes present on simulation.
        """
		n = len(self.sediment_classes())
		return n

	def mesh_data(self):
		"""
        Returns a dictionary containing information about the initial mesh data:
        how its cells and nodes are organized and related.
        """
		with open(os.path.join(self.simdir, self.mesh_path())) as mesh:
			mesh = bf(mesh, 'html.parser')
		mesh5 = self.mesh_path().replace('.ixm', '.h5')

		# Exploring Mesh .ixm to get data Adresses
		nodes = [x for x in mesh.nodes.children if x.name != None]
		to_extract = {}
		for x in nodes:
			way = x.data.text
			if x.name == 'node-ids':
				to_extract[f'nodes-{x.name}'] = way.replace('{}:'.format(mesh5), '')
			else:
				to_extract[x.name] = way.replace('{}:'.format(mesh5), '')
		cells = [x for x in mesh.cells.children if x.name != None]
		cell_data_adress = {}
		for x in cells:
			way = x.data.text
			if x.name == 'node-ids':
				to_extract[f'cells-{x.name}'] = way.replace('{}:'.format(mesh5), '')
			else:
				to_extract[x.name] = way.replace('{}:'.format(mesh5), '')

		# Getting info from Mesh .h5
		meshfile = h5py.File(os.path.join(self.simdir, mesh5), 'r')
		meshdata = {}
		for data, path in to_extract.items():
			meshdata['{}'.format(data)] = list(meshfile.get(path))
		return meshdata

	def sim_data_to_extract(self):
		"""
        Returns a dictionary contaning information needed to extract simulated results on well positions.
        'prop-data-path': dictionary relating the output file name for each timestep (KEY)
        to the path to access every output property in the corresponding  .h5 file.
        'last-data-path': path to access every property on the final timestep output file.
        'data-file': last timestep file.
        'outputpath': name of the folder containing outputfiles.
        'sed-classes': list containing sediment classes present on simulation.
        'simdir': folder containing simulation data of the "initial" project and other project inputs.
        'timesteps': list containing all simulation timesteps
        """
		with open(f'{self.simdir}/{self.output_folder()}/depouillement/data.ixm') as database:
			database = bf(database, 'html.parser')
		database = list(database.find_all('event'))
		properties_data_path = {}
		for timestep in database:
			data_file_name = timestep.attrs['name'].replace('Event', '{}'.format(self.base_file_name()))
			ts_properties = [x for x in timestep.properties.children if x.name != None]
			prop_data_path = {}
			for ts_property in ts_properties:
				ads = [x for x in ts_property.children if x.name != None]
				data_path = [x.text.replace('{}.h5:'.format(data_file_name), '') for x in ads]
				data_path = [x.replace('\n', '') for x in data_path]
				prop_data_path[ts_property['name']] = data_path

			properties_data_path[data_file_name] = prop_data_path
		last_data_path = prop_data_path
		data_file = f'{self.simdir}/{self.output_folder()}/depouillement/{data_file_name}.h5'
		sim_data_to_extract = {}
		sim_data_to_extract['prop-data-path'] = properties_data_path
		sim_data_to_extract['last-data-path'] = last_data_path
		sim_data_to_extract['data-file'] = data_file
		sim_data_to_extract['output-path'] = self.output_folder()
		sim_data_to_extract['basefilename'] = self.base_file_name()
		sim_data_to_extract['sed-classes'] = self.sediment_classes()
		sim_data_to_extract['simdir'] = self.simdir
		sim_data_to_extract['timesteps'] = self.timesteps()
		return sim_data_to_extract

class Geological_Well:
	"""
	Class used to define a Geological well object (wells used for calibration)
	"""
	
	def __init__(self, wellinfo: dict, meshcoordinates: list, markers_table):
		"""
		wellinfo: previously "table". Dictionary containing well information (name, position, facies, depth)
		meshcoordinates: list containing mesh coordinates information obtained from Dionisos_Simulation object.
		Can be obtained by 'getMeshCoordinates()' function.
		"""
		self.name = wellinfo['name']
		self.x = wellinfo['x']
		self.y = wellinfo['y']
		self.markers_table = markers_table
		self.sequence = df(data=None,columns = ['facies','depth', 'bathymetry', 'lithology'])
		self.sequence['facies'] = wellinfo['facies']
		self.sequence['depth'] = wellinfo['depth']
		if 'bathymetry' in wellinfo.keys():
			self.sequence['bathymetry'] = wellinfo['bathymetry']
		if 'lithology' in wellinfo.keys():
			self.sequence['lithology'] = wellinfo['lithology']
		#self.sequence = self.genenerateLitoBatcodes()
		self.sequence = self.applyWellMarkers()
		self.sequence = Functions.depth_to_thickness(self.sequence)
		self.sequence.dropna(subset=['facies'], inplace = True)
		self.sequence = self.dataScaling()
		self.meshcoordinates = meshcoordinates
		self.nb_data = len(self.sequence['Thickness'])
		self.facies_list = self.list_facies()
		self.geocells_matrix, self.geothick = self.well_cells() # geocells_matrix[top_position][number_of_well_data]
		self.totalthickness = self.getTotalThickness()
	
	def getTotalThickness(self):
		return self.sequence['Thickness'].sum()
	
	def list_facies(self):
		""""Returns the list of unique facies present in a well"""
		all_facies = self.sequence['facies']
		sorted_facies = sorted(list(set(all_facies)))
		return sorted_facies

	def dataScaling(self, resolution = 1.0):  # JESSICA
		"""
              Function to perform a pre-processing if a change in the geological well
              resolution is required
              Parameters:
                  thickness: List of geological well thicknesses
                  facies: List of geological well facies
                  resolution: New geological well resolution
              Returns the transformed well data (thickness and facies)
                  according to the new resolution
              """

		newThickness = []
		newSequence = self.sequence[0:0]
		if any(currentlyResolution > resolution for currentlyResolution in self.sequence['Thickness']):
			for layer in range(len(self.sequence['Thickness'])):
				if self.sequence['Thickness'][layer] / resolution > 1:
					newThickness.extend(np.ones(int(self.sequence['Thickness'][layer] // resolution)) * resolution)
					newSequence = newSequence.append([self.sequence.loc[layer]] * int(self.sequence['Thickness'][layer] // resolution))
					if self.sequence['Thickness'][layer] % resolution > 0:
						newThickness.append(round(self.sequence['Thickness'][layer] % resolution, 10))
						newSequence = newSequence.append([self.sequence.loc[layer]])
				else:
					newThickness.append(self.sequence['Thickness'][layer])
					newSequence = newSequence.append([self.sequence.loc[layer]])
			newSequence['Thickness'] = newThickness
			self.sequence = newSequence
		else:
			print("The best resolution for the geological well is already applied.")
		return (self.sequence)

	def well_cells(self):
		"""
		Method to create arrays of all the possible "cells" in a given well
		and their associated properties
		A cell is a series of data from the well with which a DionisosFlow cell 
		will be compared in the case of the stratigraphic OF
		
		: returns:
			geocells: An list of n*n lists (not exactly as it decreases as it descreases along the first n)
						of % of first facies, % of second facies...
						with n the number of well data (number of well depth)
						the first information (n) is the top position of the cell
						the second information (n) is the number of data defining the cell
			geothick: An array of size n*n (not exactly as it decreases alog the first n) with the thickness of each possible cell
						the first information (n) is the top position of the cell
						the second information (n) is the number of data defining the cell
		"""
		 
		new_table = np.asarray(pd.get_dummies(self.sequence['facies']) * np.expand_dims(self.sequence['Thickness'], axis=1))

		# initializing the matrix of all possibles cells within the well
		hiatus = np.array([0 for f in range(len(self.facies_list))])			# cells without thickness (geological unconformity)
		unique_cells = np.array([new_table[i] for i in range(self.nb_data)])	# cells containing a unique well data
		all_thicknesses = np.array(self.sequence['Thickness'])

		geothick = []
		geocells = []
		# creation of the matrix of possible cells
		for top_marker in range(self.nb_data+1):
			geothick.append([0] + [round(float(val), 2) for val in all_thicknesses[top_marker:].cumsum()])
			geocells.append(np.vstack((unique_cells[max(0, top_marker - 1)], unique_cells[top_marker:].cumsum(axis=0))))
			geocells[top_marker][1:] = np.divide(geocells[top_marker][1:].T,geothick[top_marker][1:]).T

		return geocells, geothick

	def grid_intersection(self):
		"""
        Returns a dictionary contaning:
            'relative-position': Well relative position from Simulation grid.
            'extraction-cell': nodes of the cell in the initial mesh intercepted by the well.
            This will be used in Simulated_Well object to extract the results.
        """
		origin = self.meshcoordinates[0]
		unit = self.meshcoordinates[1]
		xcoords = self.meshcoordinates[2]
		ycoords = self.meshcoordinates[3]
		nx = len(xcoords)
		ny = len(ycoords)
		[x0, y0] = origin
		x0 += (xcoords[1] - xcoords[0]) / 2
		y0 += (ycoords[1] - ycoords[0]) / 2
		origin = [x0, y0]

		position = [self.x, self.y]
		relative_position = np.array(position) - np.array(origin)

		if xcoords[-1] < relative_position[0] or ycoords[-1] < relative_position[1]:
			print(f'Problem: Well {self.name} is outside the limit of the model')
			exit()
		if relative_position[0] < 0 or relative_position[1] < 0:
			print(f'Problem: Well {self.name} is outside the limit of the model')
			exit()

		i = len([x for x in xcoords if x < relative_position[0]])
		j = len([x for x in ycoords if x < relative_position[1]])

		if i < 1 or j < 1:
			print(f'Problem: Well {self.name} is outside the grid')
			exit()

		node_ids = [(j - 1) * (nx) + i - 1, (j - 1) * (nx) + i, (j) * (nx) + i - 1, (j) * (nx) + i]
		x_fraction = (relative_position[0] - xcoords[i - 1]) / (xcoords[i] - xcoords[i - 1])
		y_fraction = (relative_position[1] - ycoords[j - 1]) / (ycoords[j] - ycoords[j - 1])
		extraction_cell = [node_ids, (x_fraction, y_fraction)]
		grid_intersection = {}
		grid_intersection['relative-position'] = relative_position
		grid_intersection['extraction-cell'] = extraction_cell
		return grid_intersection

	def applyWellMarkers(self):
		if type(self.markers_table) != str:
			zones = self.markers_table.columns.to_list()
			idx = self.markers_table.index[self.markers_table['Well'] == self.name]
			top = float(self.markers_table.loc[idx,zones[1]])
			base = float(self.markers_table.loc[idx,zones[2]])

			self.sequence.drop(self.sequence[self.sequence.depth < top].index, inplace=True)
			self.sequence.drop(self.sequence[self.sequence.depth > base].index, inplace=True)
		else:
			pass
		return self.sequence

	def genenerateLitoBatcodes(self):
		lito_refs = {'2':'4', '3':'2', '4':'3', '10':'3', '9':'2', '5':'2','7':'4','6':'3','8':'2'}
		bat_refs = {'2':'2','3':'1','4':'1','10':'3','9':'3','5':'4', '7':'1', '6':'1','8':'1'}
		self.sequence['facies'] = self.sequence.facies.astype(int)
		self.sequence['facies'] = self.sequence.facies.astype(str)
		self.sequence['lithology'] = self.sequence['facies'].replace(lito_refs)
		self.sequence['bathymetry'] =  self.sequence['facies'].replace(bat_refs)
		self.sequence['facies'] = self.sequence.facies.astype(int)
		self.sequence['lithology'] = self.sequence.lithology.astype(int)
		self.sequence['bathymetry'] = self.sequence.bathymetry.astype(int)
		self.sequence.to_excel(r'C:\\Users\\09008889978\\Desktop\\{}.xlsx'.format(self.name),index=False)
		return self.sequence

class Simulated_Well:
	"""
	Class which summarizes the information of a well so that to extract results from a DionisosFlow simulation
	"""
	def __init__(self, wellinfo:dict, relative_positions: dict, siminfo:dict, faciesdef):
		"""
		wellinfo: Dictionary containing well information (name, position, facies, depth)
		relative_positions: Well relative position from Simulation grid.
		siminfo: dictionary contaning input/output files information needed to extract simulated well results.
		faciesdef:
		"""
		self.name = wellinfo['name']
		self.x = wellinfo['x']
		self.y = wellinfo['y']
		self.relative_position = relative_positions[f'{self.name}']
		self.extraction_data = siminfo
		self.sediment_classes = siminfo['sed-classes']
		self.faciesdef = faciesdef
		self.rawdata = self.ExtractRawData()
		self.rawdata = self.facies_pre_processing()
		self.table = self.wellExtraction()
		self.nb_cells = len(self.table)
		#self.table = Functions.depth_to_thickness(self.table)
		self.table = self.reset_index()
		self.normdata = self.normalize(self.table, self.faciesdef)

	def ExtractRawData(self):
		variables = self.extraction_data['last-data-path'].keys()
		datafile = self.extraction_data['data-file']
		lastpath = self.extraction_data['last-data-path']
		propertiespath = self.extraction_data['prop-data-path']
		outputfolder = self.extraction_data['output-path']
		simdir = self.extraction_data['simdir']
		basefilename = self.extraction_data['basefilename']
		sedclasses = self.extraction_data['sed-classes']
		timesteps = self.extraction_data['timesteps']
		nsedclasses = len(sedclasses)
		cell = self.relative_position['extraction-cell']

		datatoread = h5py.File(datafile, 'r')
		welldata = {}
		for var in variables:
			varids = datatoread[lastpath[var][0]].shape[0]
			varvalues = datatoread[lastpath[var][1]].shape[0]
			ids_seq = df(datatoread[lastpath[var][0]], columns=['ids'])
			# Identifying type of variable to be read based on its length so it
			# can be read correctly.

			index = varvalues / varids

			if index == 1:
				welldata[var] = []
				print(f'Extracting results of {var}: values at each deposition age for well {self.name}.')
				for i in range(1, len(propertiespath) + 1):
					datatoread2 = h5py.File(f'{simdir}/{outputfolder}/depouillement/{basefilename}{i}.h5', 'r')
					datatoread2 = df(datatoread2.get(propertiespath[f'{basefilename}{i}'][var][1]))
					cumul = []
					for node in cell[0]:
						node = ids_seq['ids'].loc[lambda x: x == node].index[0]
						cumul.append(datatoread2.loc[node])
					result = np.array(cumul).mean()
					welldata[var].append(result)

			elif index == len(propertiespath):
				welldata[var] = []
				print(f'Extracting results of {var}: values at final age for well {self.name}.')
				varvalues = df(datatoread[lastpath[var][1]])
				cumul = []
				for node in cell[0]:
					node = ids_seq['ids'].loc[lambda x: x == node].index[0]
					cumul.append(
						np.asarray([float(varvalues.loc[node * len(propertiespath) + i]) for i in range(len(propertiespath))]))
				result = np.array(cumul).mean(axis=0)
				welldata[var] = result

			elif len(propertiespath) * nsedclasses:
				print(f'Extraction of sediment class dependent results for well {self.name}.')
				varvalues = np.array(datatoread[lastpath[var][1]])
				thickness_values = np.array(datatoread[lastpath['OutputStratigraphicThickness'][1]])
				varvalues = varvalues.reshape(int(len(varvalues) / nsedclasses), nsedclasses)
				cumul = []
				cumul_thickness = np.array([0. for x in propertiespath])
				for node in cell[0]:
					node = ids_seq['ids'].loc[lambda x: x == node].index[0]
					thickness = np.asarray(
						[thickness_values[node * len(propertiespath) + i] for i in range(len(propertiespath))])
					cumul_thickness += thickness
					concentration = np.array(varvalues[node * len(propertiespath): (node + 1) * len(propertiespath)])
					cumul.append(np.multiply(concentration.T, thickness).T)
				cumul_thickness[cumul_thickness == 0] = np.inf
				result = np.divide(np.array(cumul).sum(axis=0).T, cumul_thickness).T
				for i, litho in enumerate(sedclasses):
					welldata[f'{litho}'] = result[:, i]

		well_data_final = df(data=None, columns=welldata.keys())
		for prop, values in welldata.items():
			well_data_final[prop] = values

		# Depth "Post Processing"
		bottom = list(well_data_final['OutputBasement'])
		bottom.reverse()
		# Initializing: Only first values matter, other will be corrected below
		well_data_final['Z(Bottom)'] = [-x for x in bottom]
		well_data_final['Z(Top)'] = well_data_final['Z(Bottom)'] - well_data_final['OutputStratigraphicThickness']
		well_data_final['Z(Center)'] = (well_data_final['Z(Bottom)'] + well_data_final['Z(Top)'])/2
		for i in range(1, len(well_data_final['OutputBasement'])):
			well_data_final['Z(Bottom)'].loc[i] = well_data_final['Z(Top)'].loc[i-1]
			well_data_final['Z(Top)'].loc[i] = well_data_final['Z(Top)'].loc[i-1] - well_data_final['OutputStratigraphicThickness'].loc[i]
			well_data_final['Z(Center)'].loc[i] = (well_data_final['Z(Bottom)'].loc[i] + well_data_final['Z(Top)'].loc[i]) / 2
		well_data_final['Depth'] = -well_data_final['Z(Top)']

		# Discarding Overburden
		if len(timesteps) != len(propertiespath):
			well_data_final.drop(well_data_final.tail(1).index,inplace=True)
			print("\n\nOverburden succesfully applied. Dephts were adjusted and Overburden data discarded.")

		else:
			print("\n\nWARNING: Theres no Overburden to discard. Maybe there's a mismatch on your simulated depth data.")

		for prop, values in well_data_final.items():
			# Excluding Output from properties names
			if 'Output' in prop:
				prop_corr = prop.replace('Output', '')
				del well_data_final[prop]
			else:
				prop_corr = prop
			well_data_final[prop_corr] = values
		# Discarding Substratum
		well_data_final.drop(well_data_final.head(1).index, inplace=True)
		print(f'Substratum data (125Ma) will not be used in calibration and was successfully discarded from well {self.name}.\n')
		return well_data_final

	def normalize(self, table, faciesdef):
		 #returns the table normalized according to the table of facies definition
		for prop in faciesdef.prop_list:
			if any(table[prop]<0) or any(table[prop]>1):
				(mini, factor) = faciesdef.norms[prop]
				table[prop] = (-mini + table[prop]) / factor
		return table
		
	def distFacies(self, faciesdef):
		"""
		Functions to compute the distance between the facies of each cell with each geological facies
		and determine the facies in each cell
		Parameters:
			faciesdef: Facies_Set
		returns the table of pseudo-well data with an additional column of facies
		"""
		list_f = []
		array_dist = []
		sediment_classes = self.sediment_classes

		# in case the simulated thickness is 0, there is no sediment concentration what cause problems when computing
		# the OF. Therefore, we need to put the concentration of the first cell with thickness > 0
		# We first look in older layers and then in more recent layers
		# In case there is no thickness at all, we put the default 1 / nb of sediment classes
		total_concentration = self.normdata[sediment_classes].sum(axis=1)
		no_thickness = np.where(np.array(total_concentration) == 0.0)[0]
		thickness = np.where(np.array(total_concentration) > 0.0)[0]

		for thick in no_thickness:
			if len(thickness) == 0:
				self.normdata.iloc[thick][sediment_classes] = 1 / len(sediment_classes)
			else:
				distance = np.empty(thickness.shape)
				diff = thickness - thick
				distance[diff > 0] = diff[diff > 0]
				distance[diff < 0] = len(total_concentration) - diff[diff < 0]
				closest_time = thickness[list(distance).index(min(distance))]
				self.normdata.iloc[thick][sediment_classes] = self.normdata.iloc[closest_time][sediment_classes]

		for index, cell in self.normdata.iterrows():
			facies, distances = Functions.faciesDistances(cell, faciesdef)
			list_f.append(facies)
			array_dist.append(distances)

		self.normdata['facies'] = list_f
		self.normdata[faciesdef.faciesTab.index.values] = pd.DataFrame(np.asarray(array_dist))

	def thickness(self):
		"""
        Returns well total thickness.
        """
		thicknessdata = self.rawdata['Basement']
		thickness = thicknessdata.iloc[0] - thicknessdata.iloc[-1]
		return thickness

	def wellExtraction(self):
		add_ons = [x for x in self.sediment_classes if x not in self.faciesdef.prop_list]
		columns = add_ons + self.faciesdef.prop_list
		well = df(data=None, columns=columns)
		for prop in columns:
			well[prop] = self.rawdata[prop]
		well['depth'] = self.rawdata['Z(Bottom)']
		well['Thickness'] = self.rawdata['StratigraphicThickness']
		#Reversing df to match observed data
		well = well.iloc[::-1]
		return well

	def reset_index(self):
		self.table.reset_index(inplace=True, drop=True)
		return self.table

	def facies_pre_processing(self):
		self.rawdata['CCP'] = 0
		self.rawdata['RSandMud'] = 0
		self.rawdata['TotalCasc'] = 0
		if 'CCP' in self.faciesdef.prop_list:
			self.rawdata['CCP'] = self.rawdata['Carbo_Rud'] + self.rawdata['Carbo_Grains'] + self.rawdata['Carbo_Mud']
		if 'RSandMud' in self.faciesdef.prop_list:
			self.rawdata['RSandMud'] = (self.rawdata['Carbo_Grains'] + self.rawdata['Sand'])/(self.rawdata['Carbo_Mud']+ self.rawdata['Lutites'])
		if 'TotalCasc' in self.faciesdef.prop_list:
			self.rawdata['TotalCasc'] = self.rawdata['Gravel']+self.rawdata['Carbo_Rud']
		return self.rawdata

class Simulation_parameters:
	def __init__(self, startnumber = 0):
		"""
		Initialize the simulation parameters object.
		Parameters:
		  startnumber: Number to start the counter
		"""
		self.alloutputdata = []
		self.SCOOFoutputdata = []
		self.fullOFresults = pd.DataFrame(data=None)
		self.counter = startnumber
		self.header = []

	def updatecounter(self, increment = 1):
		"""
		Update counter number.
		Parameters:
		  increment : number to increase on the counter
		Returns the counter updated
		"""
		self.counter += increment
		return self.counter

	def saveoutput(self, replaces, OFvalue):
		"""
		Save the output on the dict.
		Parameters:
		  replaces: Uncertain parameters values to run the simulation
		  OFvalue: Objective function value for the setted parameters
		Returns the output saved
		"""
		self.header = [parameter for parameter in replaces[:,0]]
		self.header.append('OF Value')
		self.header.append('Simulation')
		candidate_solution = [value for value in replaces[:,1]]
		candidate_solution.append(OFvalue)
		candidate_solution.append(self.counter)
		self.alloutputdata.append(candidate_solution)
		return self.alloutputdata

	def saveSCOOFterms(self, SCOOF_terms):
		"""
		Save SCOOF terms on the dict.
		Parameters:
		  SCOOF_term: SCOOF value for the facies term and the thickness term
		Returns the output saved
		"""

		self.SCOOFoutputdata.append(SCOOF_terms)

		return self.SCOOFoutputdata

	def writeSCOOFtoafile(self,simdir,OF_type, well_list):
		"""
		Write the SCOOF output to a txt file.
		Parameters:
		  well_list: List of all well identification names
		Returns the txt containing all the SCOOF outputs
		"""
		results_output = f'{simdir}\\Results\\{OF_type}'
		if os.path.isdir(results_output) == False:
			os.makedirs(results_output)

		if self.counter == 0:
			self.SCOOFheader = []
			for key in well_list.keys():
				self.SCOOFheader.append(key + ' - Facies term')
				self.SCOOFheader.append(key + ' - Thickness term')
				self.SCOOFheader.append(key + ' - Beta')
			self.SCOOFheader.append('Alpha')
			self.fileDataSCOOF = PrettyTable(self.SCOOFheader)
			self.fileDataSCOOF.padding_width = 1

		self.fileDataSCOOF.add_row(self.SCOOFoutputdata[self.counter])

		with open(f'{results_output}\\SCOOF_terms_values.txt', 'w') as file:
			file.write(str(self.fileDataSCOOF))

	def writetoafile(self,simdir, OF_type, kwargs=''):
		"""
		Write the output to a txt file.
		Parameters:
		  None
		Returns the txt containing all the outputs
		"""
		results_output = f'{simdir}\\Results\\{OF_type}'
		if os.path.isdir(results_output) == False:
			os.makedirs(results_output)

		if self.counter == 0:
			self.fileData = PrettyTable(self.header)
			self.fileData.padding_width = 1
			self.alloutputdata_df = pd.DataFrame(data=[], columns=self.header)

		self.fileData.add_row(self.alloutputdata[self.counter])
		with open(f'{results_output}\\Uncertain_parameters_and_OF_values.txt', 'w') as file:
			file.write(str(self.fileData))

		round_data = pd.DataFrame(data=None,columns = self.header)
		round_data.loc[f'Sim{self.counter}',:] = self.alloutputdata[-1]
		self.alloutputdata_df = pd.concat([self.alloutputdata_df, round_data])
		if type(kwargs)!= str:
			self.alloutputdata_df.loc[f'Sim{self.counter}','Ambiguity'] = kwargs.loc['P(Ua|E)'].iloc[-1]
		with pd.ExcelWriter(f'{results_output}\\Uncertain_parameters_and_OF_values.xlsx') as writer:
			self.alloutputdata_df.to_excel(writer, sheet_name=f'OF_Values')

	def writeTotalOFresults(self, simdir, OF_type, total_OF_value, wellsOF):
		results_output = f'{simdir}\\Results\\{OF_type}'
		if os.path.isdir(results_output) == False:
			os.makedirs(results_output)
		wellsOF.index = [f'Simulation {self.counter}']
		wellsOF['TotalOF'] = total_OF_value
		self.fullOFresults = pd.concat([self.fullOFresults,wellsOF],axis=0)
		self.fullOFresults = self.fullOFresults.sort_values('TotalOF')
		with pd.ExcelWriter(f'{results_output}\\OFresults.xlsx') as writer:
				self.fullOFresults.to_excel(writer, sheet_name=f'OF_values')

	def writeTotalBeliefresults(self,simdir,OF_type,fullbelief):
		results_output = f'{simdir}\\Results\\{OF_type}\\PROOF_Stats\\Sim{self.counter}'
		if os.path.isdir(results_output) == False:
			os.makedirs(results_output)
		with pd.ExcelWriter(f'{results_output}\\Belief_WellsCombination.xlsx') as writer:
				fullbelief.to_excel(writer, sheet_name=f'Simulation_Belief')

	def writeSCOOFfacies(self, simdir, OF_type, sequences):
		results_output = f'{simdir}\\Results\\{OF_type}'
		if os.path.isdir(results_output) == False:
			os.makedirs(results_output)

		with pd.ExcelWriter(f'{results_output}\\Sim_{self.counter}_Well_Facies_Sequences.xlsx') as writer1:
			for well in sequences:
				geowell, simwell = sequences[well]
				geowell.to_excel(writer1, sheet_name=f'GeoWell{well}')
				simwell.to_excel(writer1, sheet_name=f'SimWell{well}')

	def writePROOFStats(self, simdir, OF_type, stats, evidence_combining_method):
		results_output = f'{simdir}\\Results\\{OF_type}\\PROOF_Stats\\Sim{self.counter}'
		if os.path.isdir(results_output) == False:
			os.makedirs(results_output)

		with pd.ExcelWriter(f'{results_output}\\Well_Sequences.xlsx') as writer1:
			for well in stats:
				if evidence_combining_method == 'belief':
					geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends, df_prob_f], [T_Values,F_Values], full_well_belief, belief = 	stats[well]
				elif evidence_combining_method == 'mean':
					[geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends, df_prob_f], [T_Values, F_Values]] = stats[well]
				geowell.loc['ID',:] = 'GeoWell'
				simwell.loc['ID',:] = 'SimWell'
				welltable = pd.concat([geowell,simwell],axis=1)
				welltable.to_excel(writer1, sheet_name=f'{well}')

		with pd.ExcelWriter(f'{results_output}\\Well_Stats.xlsx') as writer2:
			for well in stats:
				if evidence_combining_method == 'belief':
					geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends, df_prob_f], [T_Values, F_Values], full_well_belief, belief = stats[well]
				elif evidence_combining_method == 'mean':
					[geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends, df_prob_f], [T_Values, F_Values]] = stats[well]
				geowellstats.loc['ID', :] = 'GeoWell'
				geowellstats.loc['-', :] = '-'
				simwellstats.loc['ID', :] = 'SimWell'
				simwellstats.loc['-', :] = '-'
				T_F_values = pd.concat([T_Values, F_Values], axis=0)
				T_F_values.loc['-', :] = '-'
				df_prob_t.index = ['To_Be_Prob (Mean)']
				df_prob_f.index = ['To_Be_Prob_(Var)']
				df_prob_t_tends.index = ['To_Be_Prob_(Tends Mean)']
				probabilities = pd.concat([df_prob_t, df_prob_t_tends, df_prob_f], axis=0)
				probabilities.loc['-', :] = '-'
				probabilities.fillna('-', inplace=True)
				statstable = pd.concat([geowellstats, simwellstats, T_F_values, probabilities], axis=0)
				statstable.to_excel(writer2, sheet_name=f'{well}')

		if evidence_combining_method == 'belief':
			with pd.ExcelWriter(f'{results_output}\\Well_evidences_belief.xlsx') as writer4:
				for well in stats:
					geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends,  df_prob_f], [T_Values, F_Values], full_well_belief, belief = stats[well]
					full_well_belief.to_excel(writer4, sheet_name=f'{well}')
