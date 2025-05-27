#Import libraries
import Code_CW.Functions as Functions
import pandas as pd
from pandas import DataFrame as df
from bs4 import BeautifulSoup as bf
import os
import pathlib
import lasio
import h5py
import numpy as np
from prettytable import PrettyTable
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt
import matplotlib
import warnings
import pathvalidate

class FaciesAnalysis:
	def __init__(self, root, wells_dir, markers_file, extract_columns, domain, color_reference_file):
		warnings.simplefilter(action='ignore', category=FutureWarning)
		self.root = root
		self.domain = domain
		self.setdomain()
		self.wells_dir = wells_dir
		self.well_files = Functions.getWellFiles(self.wells_dir)
		self.markers_file = markers_file
		self.markers_table = self.readMarkersFile(markers_file)
		self.faciescolor = self.set_color_map(color_reference_file)
		self.extract_columns = extract_columns
		self.welldata = self.getGeoWellData(self.well_files, self. markers_table, self.extract_columns)
		self.facies = self.registerAllFacies(self.welldata)
		self.facieslist = [f'facies {f}' for f in self.facies]
		self.facieslist.append('total')
		self.wellFaciesAccounting(self.facies)
		self.faciesinfo = self.FaciesInfo(self.welldata,self.facies)
		self.generalinfo = self.GeneralInfo(self.welldata, self.markers_table )
		self.totalinfo, self.wellscatterdata = self.TotalInfo(self.welldata, self.faciesinfo, self.generalinfo)
		self.registerXLSX(self.root, self.totalinfo)

	def readMarkersFile(self, markers_file):
		if markers_file != '':
			f = open(markers_file, 'r')
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
				print('Error: the file is not a file of definition of well/zone markers')
			else:
				try:
					nb_wells = int(f.readline())
				except Exception as e:
					print('Exception:', e)

			f.close()

			markers_table = pd.read_csv(markers_file, delim_whitespace=True, header=None, skiprows=(nb_zones + 4))
			markers_table.columns = zones
		else:
			markers_table = ''
		return markers_table

	def set_color_map(self, color_reference_file):
		if color_reference_file != '':
			colors_ref = pd.read_csv(color_reference_file, delim_whitespace=True, header=0)
			colors_ref.sort_values(by='Facies',inplace =True)
			for i in range(0, len(colors_ref.Facies)):
				colors_ref.Facies.loc[i] = 'facies '+ str(colors_ref.Facies.loc[i])
			colors_ref.set_index('Facies', inplace=True)
			colors_ref.dropna(inplace=True)
		else:
			colors_ref = ''
		return colors_ref

	def readWellFile(self, well, extract_columns):
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
		# Try to find well properties on .las file. If not found, returns a ERROR message.
		# WELL NAME
		try:
			wellinfo['name'] = file.well.WELL.value
		except AttributeError:
			print("!!ERROR: Could not find well name on .las file.!!")
		# WELL COORDINATES
		# Checking if there's _dev files on well list:
		check_file = well.replace('.las', '_dev')
		if os.path.exists(check_file) == True:
			try:
				wellinfo['x'] = Functions.getPropFromTxtFile(check_file, 'X-COORDINATE')
				wellinfo['y'] = Functions.getPropFromTxtFile(check_file, 'Y-COORDINATE')
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
					print(
						f'!!ERROR: Theres no X,Y coordinates information on .las file.\n Check the {wellinfo["name"]} .las file.')
					exit()
		# WELL DEPTH
		try:
			wellinfo['depth'] = file.curves[f'{depth}'].data
		except KeyError:
			print(
				f'!!ERROR: Theres no column named {depth}.\n Check the DEPTH column name to be extracted on GeoWell {wellinfo["name"]}!!')
		# WELL FACIES
		try:
			wellinfo['facies'] = file.curves[f'{facies}'].data
		except KeyError:
			print(
				f'!!ERROR: Theres no column named {facies}.\n Check the FACIES column name to be extracted on GeoWell {wellinfo["name"]}!!')
		# WELL BATHYMETRY
		try:
			wellinfo['bathymetry'] = file.curves[f'{bathymetry}'].data
		except KeyError:
			print(
				f'!!ERROR: Theres no column named {bathymetry}.\n Check the BATHYMETRY column name to be extracted on GeoWell {wellinfo["name"]}!!')
		# WELL LITHOLOGY
		try:
			wellinfo['lithology'] = file.curves[f'{lithology}'].data
		except KeyError:
			print(
				f'Theres no column named {lithology}. \n Check the LITHOLOGY column name to be extracted on GeoWell {wellinfo["name"]}')
		return wellinfo

	def getGeoWellData(self, files, markers_table, extract_columns):
		welldata = {}
		for well in files:
			wellinfo = Functions.readWellFile(well, extract_columns)
			welldata[f"{wellinfo['name']}"] = FA_GeoWell(wellinfo, markers_table)
		return welldata

	def registerAllFacies(self, welldata):
		allfacies = []
		for well,data in welldata.items():
			data.facies.dropna(axis=0, inplace=True)
			for facies in data.facies['facies']:
				if facies not in allfacies:
					allfacies.append(facies)
				else:
					pass
		allfacies = np.sort(allfacies)
		allfacies = list(allfacies)
		return allfacies

	def wellFaciesAccounting(self, facies):
		for well,data in self.welldata.items():
			faciesthickness = df(data=None, index = ['ZOI'], columns = [self.facieslist])
			for f in facies:
				#Well
				zoitemp = data.sequence.loc[data.sequence['facies'] == f].sum()
				zoitemp = zoitemp['Thickness']
				faciesthickness.iloc[0].loc[f'facies {f}'] = zoitemp
			faciesthickness.iloc[0].loc['total'] = np.sum(faciesthickness.loc['ZOI'][:])
			faciesthickness.fillna(0,inplace=True)
			data.faciesthickness = faciesthickness

	def FaciesInfo(self, welldata, facies):
		faciesinfo = {}
		for well, data in welldata.items():
			faciesproportion = df(data=None, index=['ZOI'], columns=[self.facieslist])
			for f in facies:
				faciesproportion.loc['ZOI'][f'facies {f}'] = float(data.faciesthickness.iloc[0].loc[f'facies {f}']) / float(data.faciesthickness.iloc[0].loc['total'])
			faciesproportion.fillna(0, inplace=True)
			faciesproportion.loc['ZOI','total'] = np.sum(faciesproportion.loc['ZOI'][:])
			wellfaciesinfo = {'Thickness': data.faciesthickness, 'Proportion': faciesproportion}
			faciesinfo[well] = wellfaciesinfo
		return faciesinfo

	def setdomain(self):
		self.xo = self.domain[0]
		self.xlen = self.domain[1]
		self.yo = self.domain[2]
		self.ylen = self.domain[3]

	def defineWellLocation(self, welldata, well):
		wellx = welldata[well].x
		welly = welldata[well].y
		x_mid_axis = self.xo+(self.xlen/2)
		y_mid_axis = self.yo + (self.ylen / 2)
		if wellx == x_mid_axis and welly == y_mid_axis:
			welllocation = 'Center'
		elif wellx == x_mid_axis and welly != y_mid_axis:
			if welly < y_mid_axis:
				welllocation = 'South'
			elif welly > y_mid_axis:
				welllocation = 'North'
		elif welly == y_mid_axis and wellx != x_mid_axis:
			if wellx < x_mid_axis:
				welllocation = 'West'
			elif wellx > x_mid_axis:
				welllocation = 'East'
		elif wellx < x_mid_axis and welly < y_mid_axis:
			welllocation = 'SouthWest'
		elif wellx < x_mid_axis and welly > y_mid_axis:
			welllocation = 'NorthWest'
		elif wellx > x_mid_axis and welly < y_mid_axis:
			welllocation = 'SouthEast'
		elif wellx > x_mid_axis and welly > y_mid_axis:
			welllocation = 'NorthEast'
		return welllocation

	def GeneralInfo(self, welldata, markers_table):
		generalinfo = {}
		for well,data in welldata.items():
			welllocation = self.defineWellLocation(welldata, well)
			wellinfo = df(np.zeros((4, 1)), index=['Location', 'East Coordinate', 'North Coordinate', 'Well Thickness'],
						  columns=['Well'])
			wellinfo.loc['Location'] = welllocation
			wellinfo.loc['East Coordinate'] = data.x
			wellinfo.loc['North Coordinate'] = data.y
			wellinfo.loc['Well Thickness'] = np.sum(data.sequence.Thickness[:])
			markers_temp = markers_table.copy(deep=True)
			markers_temp.drop('Well', inplace=True, axis=1)
			idx = self.markers_table.index[self.markers_table['Well'] == int(well)]
			for marker in markers_temp:
					MD = float(markers_table.loc[idx,marker])
					wellinfo.loc[f'{marker} MD'] = MD
			generalinfo[well] = wellinfo
		return generalinfo

	def TotalInfo(self, welldata, faciesinfo, generalinfo):
		wellscatterdata = df(data=None, index=[f'Well {x}' for x in self.welldata.keys()],
							 columns=['East', 'North', 'Dominant Facies', 'Size', 'Facies'])
		totalinfo = {}
		for well,data in welldata.items():
			totalinfo[f'Well {well}'] = {'WELL': generalinfo[well], 'FACIES': faciesinfo[well]}
			wellscatterdata.loc['Well {}'.format(well)]['East'] = generalinfo[well]['Well']['East Coordinate']
			wellscatterdata.loc['Well {}'.format(well)]['North'] = generalinfo[well]['Well']['North Coordinate']
		return totalinfo, wellscatterdata

	def registerXLSX(self, root, totalinfo):
		if os.path.isdir(f'{root}\Output') == False:
			os.makedirs(f'{root}\Output')

		self.allfacies_to_plot = self.facieslist + ['X', 'Y']

		facies_thickness_to_excel = df(data=None, index = totalinfo.keys(), columns = [self.allfacies_to_plot])
		for well, data in totalinfo.items():
			facies_thickness_to_excel.loc[well]['X'] =data['WELL']['Well']['East Coordinate']
			facies_thickness_to_excel.loc[well]['Y'] = data['WELL']['Well']['North Coordinate']
			for f in self.facieslist:
				facies_thickness_to_excel.loc[well][f] = data['FACIES']['Thickness'].loc['ZOI'][f'{f}']

		facies_proportions_to_excel = df(data=None, index = totalinfo.keys(), columns = [self.allfacies_to_plot])
		for well, data in totalinfo.items():
			facies_proportions_to_excel.loc[well]['X'] =data['WELL']['Well']['East Coordinate']
			facies_proportions_to_excel.loc[well]['Y'] = data['WELL']['Well']['North Coordinate']
			for f in self.facieslist:
				facies_proportions_to_excel.loc[well][f] = data['FACIES']['Proportion'].loc['ZOI'][f'{f}']

		with pd.ExcelWriter(f'{root}\Output\\WellFaciesValues.xlsx') as writer:
			facies_thickness_to_excel.to_excel(writer,sheet_name='Thickness Data')
			facies_proportions_to_excel.to_excel(writer,sheet_name='Proportions Data')
		return

	def wellPieCharts(self, totalinfo, faciescolor):
		if os.path.isdir(f'{self.root}\\Output\\WellPieCharts') == False:
			os.makedirs(f'{self.root}\\Output\\WellPieCharts')

		for well, data in totalinfo.items():
			data_to_plot = data['FACIES']['Proportion']
			data_to_plot.drop('total', axis = 1, inplace=True)

			for x in data_to_plot.columns:
				if data_to_plot.iloc[0].loc[x] == 0:
					data_to_plot.drop(x, axis=1, inplace = True)

			wellcolors = [faciescolor.loc[x] for x in faciescolor.index if x in data_to_plot.columns]
			labels = [x for x in data_to_plot.columns]
			# gs = gridspec.GridSpec(1,3, width_ratios=[1,1,1], height_ratios=[1])
			# ax1 = plt.subplot(gs[0,0])

			def formatpctvalues(data):
				pct = '{:.1f}%'.format(data)
				return pct

			plt.pie(data_to_plot.iloc[0], labels=labels, colors=wellcolors, autopct= formatpctvalues, pctdistance=0.65, startangle=0,
					textprops=dict(size=8,fontweight = 'bold'))

			plt.title(f'{well} Facies Proportions')
			plt.savefig(f'{self.root}\\Output\\WellPieCharts\\{well}', transparent=True)
			plt.clf()

	def getImage(self,path):
		return OffsetImage(plt.imread(path), zoom=0.12)

	def scatterpie(self, totalinfo):
		path = f'{self.root}\\Output\\WellPieCharts'
		wellcoord = df(data=None, columns=['X', 'Y', 'PATH'])
		pies =[]
		for r, d, f in os.walk(path):
			for file in f:
					pies.append(os.path.join(r, file))

		for pie in pies:
			base = os.path.basename(pie)
			wellpie = os.path.splitext(base)
			wellpie = f'{wellpie[0]}'

			wellcoord.loc[f'{format(wellpie)}', 'X'] = float(totalinfo[wellpie]['WELL'].loc['East Coordinate'])
			wellcoord.loc[f'{format(wellpie)}', 'Y'] = float(totalinfo[wellpie]['WELL'].loc['North Coordinate'])
			wellcoord.loc[f'{format(wellpie)}', 'PATH'] = pie
		fig, ax = plt.subplots(dpi=1200)
		ax.scatter(wellcoord.X, wellcoord.Y, s=0.02)
		for i in wellcoord.index:
			ab = AnnotationBbox(self.getImage(wellcoord.PATH.loc[i]), (wellcoord.X.loc[i], wellcoord.Y.loc[i]),
								frameon=False, pad=1.3)
			ax.add_artist(ab)
		plt.xlabel('East')
		plt.xlim(self.xo, self.xo+self.xlen)
		plt.ylabel('North')
		plt.ylim(self.yo, self.yo+self.ylen)
		plt.savefig(f'{self.root}\\Output\\ScatterPie.png')
		plt.clf()

class FA_GeoWell:
	def __init__(self, wellinfo, markers_table):
		self.name = wellinfo['name']
		self.x = wellinfo['x']
		self.y = wellinfo['y']
		self.rawdata = {'facies': wellinfo['facies'], 'depth': wellinfo['depth'], 'lithology':wellinfo['lithology']}
		self.markers_table = markers_table
		self.sequence = df(self.rawdata)
		self.sequence = Functions.depth_to_thickness(self.sequence)
		self.sequence = self.applyWellMarkers(self.markers_table)
		self.facies = self.getWellFacies(self.sequence)
		self.thickness = np.sum(self.sequence['Thickness'])

	def applyWellMarkers(self, markers_table):
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

	def getWellFacies(self, sequence):
		faciesid = sequence.drop_duplicates('facies')
		faciesid = df((faciesid.facies))
		faciesid.reset_index(drop=True,inplace=True)
		return faciesid

#Import libraries
import Code_CW.Functions as Functions
import pandas as pd
from pandas import DataFrame as df
from bs4 import BeautifulSoup as bf
import os
import pathlib
import lasio
import h5py
import numpy as np
from prettytable import PrettyTable
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt
import matplotlib
import warnings

class SensitivityAnalysis:
	def __init__(self,case, OF, simdir, origin, welldir, extract_columns, facies_definition, simulation, enginefolder,
				 alpha = 0, beta = 0, internal_optimization_method='', litho_class_ids = '', carbo_ids = '', mud_ids='', sand_ids='', gravel_ids ='',
				 markers_table = '',reference_value = [], expdesign = '', uncertain_params = '', evidence_combining_method = 'mean', total_OF_computing_method = 'mean'):
		self.case = case
		self.OF = OF
		self.OF_type = OF.OF_type
		self.simdir = simdir
		self.origin = origin
		self.extract_columns = extract_columns
		self.facies_definition = facies_definition
		self.simulation = simulation
		self.enginefolder = enginefolder
		self.OF.reference_value = reference_value
		self.markers_table = markers_table
		self.set_alpha = alpha
		self.set_beta = beta
		self.litho_class_ids = litho_class_ids
		self.carbo_ids = carbo_ids
		self.mud_ids = mud_ids
		self.sand_ids = sand_ids
		self.gravel_ids = gravel_ids
		self.expdesign = expdesign
		self.uncertain_params = uncertain_params
		self.scenariolist = self.getSimList()
		self.well_files = Functions.getWellFiles(welldir)
		self.internal_optimization_method = internal_optimization_method
		self.evidence_combining_method = evidence_combining_method
		self.total_OF_method = total_OF_computing_method
		pass

	def computeSensitivityAnalysis(self):
		scenariosOF = {}
		if self.OF_type =='P':
			self.setReferenceSolution()
		for scenario in self.scenariolist:
			self.extractScenarioResults(scenario)
			self.OF.alpha = self.set_alpha if self.set_alpha != 0 else (1 / np.sqrt(self.facies_definition.nb_prop))
			self.OF.SCOOF_terms = []
			scenario_name = f'__Simu {scenario.split("_")[-1]}__'
			if self.OF_type == 'P':
				if self.total_OF_method == 'belief':
					total_OF_value, wells_OF, fullstats, total_OF_belief = self.computeScenarioOF()
					scenariosOF[scenario_name] = [total_OF_value, wells_OF, fullstats, total_OF_belief]
				elif self.total_OF_method == 'mean':
					total_OF_value, wells_OF, fullstats = self.computeScenarioOF()
					scenariosOF[scenario_name] = [total_OF_value, wells_OF, fullstats]
			else:
				total_OF_value, wells_OF, fullstats, total_OF_belief = self.computeScenarioOF()
				scenariosOF[scenario_name] = [total_OF_value, wells_OF, fullstats, total_OF_belief]
			print(f'{scenario_name} OF value: {total_OF_value}\n\n')
		return scenariosOF

	def computeScenarioOF(self):
		wells_belief = pd.DataFrame(np.zeros((3, len(self.OF.geowell.keys()))), columns=self.OF.geowell.keys())
		wells_OF = pd.DataFrame(np.zeros((1, len(self.OF.geowell.keys()))), columns=self.OF.geowell.keys(), index = ['Well OF Value'])
		wellstats = {}
		#COMPUTE EVERY SCENARIO FO AND REGISTER ITS VALUE, WELL BY WELL
		for key in self.OF.geowell.keys():
			self.OF.simwell = self.OF.simwell_data[key]
			self.OF.well = self.OF.geowell.get(key)
			if self.OF_type == 'S':
				self.OF.beta = self.set_beta if self.set_beta != 0 else (1 / self.OF.well.totalthickness)
				self.OF.simwell.distFacies(self.facies_definition)
				well_OF = self.OF.computeWellOF()
				sequences = [self.OF.well.sequence, self.OF.simwell.table]
				wellstats[key] = sequences

			if self.OF_type == 'P':
				if self.OF.reference_value != []:
					self.OF.well = self.ref_solution.get(key)
				well_OF, wellstats[key] = self.OF.computeWellOF_fromsamples()
				if self.evidence_combining_method =='belief':
					wells_belief[key] = wellstats[key][-1]
			wells_OF[key] = well_OF

		if self.OF_type == 'S':
			total_OF_value = np.mean(wells_OF.iloc[0, :])
			print("Partial OF value: ", total_OF_value)
			self.OF.SCOOF_terms.append(self.OF.alpha)
			self.OF.simulation.saveSCOOFterms(self.OF.SCOOF_terms)
			self.OF.simulation.writeSCOOFtoafile(self.simdir,self.OF_type, self.OF.geowell)
			self.OF.simulation.updatecounter()
			return total_OF_value, wells_OF, wellstats, 0

		elif self.OF_type == 'P':
			if (self.total_OF_method =='belief' and self.evidence_combining_method =='mean'):
				wells_OF_to_belief = wells_OF.copy(deep=True)
				wells_OF_to_belief.loc['Well OF Value', :] = 1 - wells_OF_to_belief.loc['Well OF Value', :]
				total_OF_value, probabilities, full_well_belief = self.OF.BeliefFunctionEvidenceComputing(wells_OF_to_belief)
				print("Partial OF value: ", total_OF_value)
				return total_OF_value, wells_OF, wellstats, full_well_belief
			#wells_OF_to_belief = wells_OF.copy(deep=True)
			#wells_OF_to_belief = 1 - wells_OF_to_belief.loc['Well OF Value', :]
			elif (self.total_OF_method == 'belief' and self.evidence_combining_method =='belief'):
				total_OF_value, total_OF_belief = self.OF.BeliefFunctionTotalOFComputing(wells_belief)
				print("Partial OF value: ", total_OF_value)
				return total_OF_value, wells_OF, wellstats, total_OF_belief
			elif (self.total_OF_method == 'mean' and self.evidence_combining_method =='mean'):
				total_OF_value = np.mean(wells_OF.loc['Well OF Value'])
				print("Partial OF value: ", total_OF_value)
				return total_OF_value, wells_OF, wellstats

	def extractScenarioResults(self,scenario):
		#EXCTRACTS EACH SCENARIO RESULT
		scenariodir = os.path.join(self.simdir, 'Simulations', scenario)
		arcfilebase = Functions.getArcInput(scenariodir)
		self.OF.arcinputfile = arcfilebase
		self.OF.set_objective_function(self.case, scenariodir, arcfilebase, self.origin, self.well_files, self.extract_columns, self.facies_definition, self.simulation,
									   markers_table = self.markers_table,
									   litho_class_ids = self.litho_class_ids, carbo_ids=self.carbo_ids, mud_ids = self.mud_ids, sand_ids= self.sand_ids, gravel_ids=self.gravel_ids,
									   internal_optimization_method = self.internal_optimization_method,
									   evidence_combining_method=self.evidence_combining_method, total_OF_computing_method=self.total_OF_method)
		if self.OF_type == 'P':
			self.OF.set_lithologies()
		self.OF.extract_simulation_data()
		self.OF.extract_obs_well_data()
		self.OF.extract_sim_well_data()
		return self.OF

	def getSimList(self):
		#GETS LIST OF CENARIOS FROM SIMULATION DIRECTORY
		folders = next(os.walk(os.path.join(self.simdir, 'Simulations')))
		scenariolist = folders[1]
		return scenariolist

	def setReferenceSolution(self):
		#WHEN WORKING WITH MIXED CASE AND USING PROOF, IT IS NECESSARY TO ESTABLISH A REFERENCE MODEL
		if self.OF.reference_value != [] and self.simulation.counter == 0:
			arcfilebase = Functions.getArcInput(self.simdir)
			ref_replace = self.OF.reference_value
			ref_arcinput = Functions.replaceInversibleParameters(ref_replace, self.simdir, arcfilebase)
			self.OF.arcinputfile = ref_arcinput
			Functions.runDionisosSim(self.simdir, self.enginefolder, ref_arcinput)
			self.OF.set_objective_function(self.case, self.simdir, arcfilebase, self.origin, self.well_files,
										   self.extract_columns, self.facies_definition, self.simulation)
			self.OF.extract_simulation_data()
			self.OF.extract_obs_well_data()
			self.OF.ref_solution = Functions.getSimWellData(self.well_files, self.simdir,
													Functions.getGeoWellRelativePosition(self.OF.geowell),
													Functions.toWellExtraction(self.simdir, ref_arcinput, self.origin,
																			   self.extract_columns), self.facies_definition)
			return
		else:
			pass

	def writeSensResults(self, OF_values):
		print('Registering Sensitivity Analysis results.')
		#Create Results Folder
		if os.path.isdir(f'{self.simdir}\Sens_Results') == False:
			os.makedirs(f'{self.simdir}\Sens_Results')
		#Create PROOF Results Folder
		if self.OF_type =='P':
			print(f'Minimum OF value obtained in {min(OF_values, key=OF_values.get)} through Probabilistic Objective Function (PROOF) methods. \nOF value = {min(OF_values.values())[0]}')
			oFfolder = 'PROOF'
			if os.path.isdir(f'{self.simdir}\Sens_Results\{oFfolder}') == False:
				os.makedirs(f'{self.simdir}\Sens_Results\{oFfolder}')
		#Create SCOOF Results Folder
		elif self.OF_type =='S':
			oFfolder = 'SCOOF'
			print(f'Minimum OF value obtained in {min(OF_values, key=OF_values.get)} through Stratigraphic Correlation Objective Function (SCOOF) methods. \nOF value = {min(OF_values.values())[0]}')
			if os.path.isdir(f'{self.simdir}\Sens_Results\{oFfolder}') == False:
				os.makedirs(f'{self.simdir}\Sens_Results\{oFfolder}')

		if self.expdesign != '':
			self.experimental = self.readExpDesign()
			uncertainty = self.readUncertaintyranges()
			uncertainty = self.uncertainty_values_processing(uncertainty)
			params_to_plot = self.parameters_adjust_to_plot(uncertainty)
			param_table, param_table_to_plot = self.getparamvalues(uncertainty, params_to_plot)
		else:
			param_table = []
			param_table_to_plot = []

		ofs_to_excel, wellofs_results = self.writeOFsTable(OF_values,oFfolder, param_table, param_table_to_plot)
		self.writeTXTFile(oFfolder, OF_values)
		self.plotGlobalSimsOFs(ofs_to_excel, oFfolder, param_table_to_plot, wellofs_results)
		self.plotWellSimsOFs(oFfolder, param_table_to_plot, wellofs_results)
		if self.OF_type == 'P':
			self.writePROOFstats(oFfolder, OF_values, self.evidence_combining_method)
			if self.total_OF_method =='belief':
				self.writeTotalBeliefresults(oFfolder, OF_values)
		else:
			if self.expdesign != '':
				facies_terms_to_plot, thickness_terms_to_plot = self.processingSCOOFcontribution()
				self.plotGlobalContributions(oFfolder, ofs_to_excel, param_table_to_plot, facies_terms_to_plot, thickness_terms_to_plot)
				self.plotWellsContributions(oFfolder, ofs_to_excel, param_table_to_plot, facies_terms_to_plot, thickness_terms_to_plot)

	def writeTotalBeliefresults(self,oFfolder,OF_values):
		if os.path.isdir(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats') == False:
			os.makedirs(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats')
		for simu, result in OF_values.items():
			fullbelief = result[-1]
			with pd.ExcelWriter(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats\\Simu_{simu}\\Belief_WellsCombination.xlsx') as writer:
					fullbelief.to_excel(writer, sheet_name=f'Simulation_Belief')

	def readUncertaintyranges(self):
		#READS UNCERTAIN PARAMETERS FILE
		uncertainty = pd.read_csv(self.uncertain_params, encoding= 'unicode_escape', sep = '\t')
		return uncertainty

	def readExpDesign(self):
		# READS EXPERIMENTAL DESIGN FILE
		experimental = pd.read_excel(self.expdesign)
		experimental.drop('Simulation', axis=1, inplace=True)
		experimental.drop('Third party results', axis=1, inplace=True)
		experimental.drop([x for x in experimental.columns if ('Design' or 'Exp') in x], axis=1, inplace=True)
		return experimental

	def plotGlobalSimsOFs(self, ofs_to_excel, oFfolder, param_table_to_plot, wellofs_results):
		for param in param_table_to_plot:
			if os.path.isdir(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}') == False:
				os.makedirs(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}')
			ofs_to_excel = ofs_to_excel.sort_values(by=f'{param}')
			sims = []
			for sim in ofs_to_excel['simulations']:
				emp_str = ""
				for m in sim:
					if m.isdigit():
						emp_str = emp_str + m
				sims.append(int(emp_str))

			fig, ax = plt.subplots()
			fig.set_size_inches(12, 5)
			plt.plot(ofs_to_excel[param], ofs_to_excel['OF value'], marker='o', linewidth=1.5, markerfacecolor='k',
					 color='k')
			plt.xlabel(f'{param}')
			plt.ylabel('OF value')
			plt.title(f'OF_values versus {param}')
			x = ofs_to_excel[param].reset_index(drop=True)
			y = ofs_to_excel['OF value'].reset_index(drop=True)
			n = [str(t) for t in sims]
			for i, txt in enumerate(n):
				ax.annotate(txt, (x[i], y[i]), horizontalalignment='center', verticalalignment='top', xytext=(0, 15),
							size=10, textcoords='offset points')

			plt.savefig(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}\\{oFfolder}_values.png')
			plt.close()

	def uncertainty_values_processing(self,uncertainty):
		#PROCESSES UNCERTAINTY INFOS SO IT CAN BE PROPERLY USED AND REGISTERED
		i = 0
		n = len(uncertainty)
		while i < n:
			if uncertainty.loc[i].loc['Type'] == 'Real':
				if float(uncertainty.loc[i].loc['DEF_VALUES']) == 0 and float(uncertainty.loc[i].loc['MIN']) == -1 and float(uncertainty.loc[i].loc['MAX']) == 1:
					idx = i
					i += 1
					while uncertainty.loc[i].loc['Variables'][0:2] =='S+' or uncertainty.loc[i].loc['Variables'][0:2] =='M+':
						if uncertainty.loc[i].loc['Variables'][0] =='S':
							uncertainty = uncertainty.drop(i, axis=0)
						elif uncertainty.loc[i].loc['Variables'][0] =='M':
							uncertainty.loc[idx,'DEF_VALUES'] = float(uncertainty.loc[i,'DEF_VALUES'])
							uncertainty.loc[idx,'MIN'] = float(uncertainty.loc[i,'MIN'])
							uncertainty.loc[idx,'MAX'] = float(uncertainty.loc[i,'MAX'])
							uncertainty = uncertainty.drop(i, axis=0)
						i += 1
						if i >= n:
							break
				else:
					uncertainty.loc[i].loc['DEF_VALUES'] = float(uncertainty.loc[i].loc['DEF_VALUES'])
					uncertainty.loc[i].loc['MIN'] = float(uncertainty.loc[i].loc['MIN'])
					uncertainty.loc[i].loc['MAX'] = float(uncertainty.loc[i].loc['MAX'])
					i+= 1
			elif uncertainty.loc[i].loc['Type'] == 'Curve':
				ts_values = uncertainty.loc[i].loc['DEF_VALUES'].split('],[')
				ts_values.pop(0)
				ts_values = ts_values[0].replace(']]', '').split(',')
				uncertainty.loc[i,'DEF_VALUES'] = ','.join(ts_values)
				ts_values = uncertainty.loc[i].loc['MIN'].split('],[')
				ts_values.pop(0)
				ts_values = ts_values[0].replace(']]', '').split(',')
				uncertainty.loc[i,'MIN'] = ','.join(ts_values)
				ts_values = uncertainty.loc[i].loc['MAX'].split('],[')
				ts_values.pop(0)
				ts_values = ts_values[0].replace(']]', '').split(',')
				uncertainty.loc[i,'MAX'] = ','.join(ts_values)
				i += 1

			elif uncertainty.loc[i]['Type'] =='Map':
				ts_values = uncertainty.loc[i].loc['DEF_VALUES'].split("],[")
				ts_values = ts_values[0].replace('[[','')
				initial_value = float(ts_values.split(',')[0])
				uncertainty.at[i,'DEF_VALUES'] = 1
				ts_values = uncertainty.loc[i].loc['MIN'].split("],[")
				ts_values = ts_values[0].replace('[[', '')
				minimum = float(ts_values.split(',')[0])
				uncertainty.at[i,'MIN'] = minimum/initial_value
				ts_values = uncertainty.loc[i].loc['MAX'].split("],[")
				ts_values = ts_values[0].replace('[[', '')
				maximum = float(ts_values.split(',')[0])
				uncertainty.at[i,'MAX'] = maximum/initial_value
				i += 1
			else:
				i += 1
		uncertainty.reset_index(inplace=True,drop=True)
		return uncertainty

	def parameters_adjust_to_plot(self, uncertainty):
		plot_params = uncertainty.copy(deep=True)
		plot_params = plot_params.drop('Unit',axis=1)
		for i in range(0, len(plot_params)):
			if plot_params.loc[i].loc['Type'] =='Map' or plot_params.loc[i].loc['Type'] =='Real':
				continue
			else:
				inits = plot_params.loc[i].loc['DEF_VALUES'].split(',')
				mins = plot_params.loc[i].loc['MIN'].split(',')
				maxs = plot_params.loc[i].loc['MAX'].split(',')
				initial_value = float(inits[-1])
				if plot_params.loc[i].loc['Type'] == 'Curve':
					if float(mins[-1]) != initial_value:
						plot_params.loc[i,'DEF_VALUES'] = 1
						minimum = float(mins[-1])
						maximum = float(maxs[-1])
					elif float(mins[0]) != initial_value and float(mins[0]) != 0:
						initial_value = float(inits[0])
						plot_params.loc[i,'DEF_VALUES'] = 1
						minimum = float(mins[0])
						maximum = float(maxs[0])
					elif float(mins[0]) == initial_value and float(mins[-1]) == 0:
						initial_value = float(inits[0])
						plot_params.loc[i, 'DEF_VALUES'] = 1
						minimum = 0
						maximum = float(maxs[0])
					plot_params.loc[i,'MIN'] = float(minimum / initial_value)
					plot_params.loc[i,'MAX'] = float(maximum / initial_value)

		plot_params.reset_index(inplace=True, drop=True)
		return plot_params

	def params_def(self,imin, init, imax):
		if type(init) == list:
			params = df(data=None, columns=['min', 'init', 'max'])
			params['min'] = imin
			params['init'] = init
			params['max'] = imax
		else:
			params = df(data=None, index=[0])
			params['min'] = imin
			params['init'] = init
			params['max'] = imax
		return params

	def interpol(self,params, index):
		x = []
		for i in range(len(params)):
			minimo = float(params['min'][i])
			initial = float(params['init'][i])
			maximo = float(params['max'][i])
			if index > 0:
				value = 1
				y = -((maximo - initial) * (value - index) / value) + maximo
				x.append(y)
			elif index < 0:
				value = -1
				y = -(((initial - minimo) * (-index)) / (-value)) + initial
				x.append(y)
			elif index == 0:
				x.append(initial)
		return x

	def getparamvalues(self, uncertainty, plot_params):
		paramvalues = pd.DataFrame(np.zeros([len(self.experimental),len(self.experimental.columns)]),columns = self.experimental.columns)
		plot_paramvalues = pd.DataFrame(np.zeros([len(self.experimental),len(self.experimental.columns)]),columns = self.experimental.columns)
		for param in self.experimental.columns:
			full_param_value = []
			full_param_value_to_plot = []
			param_aux = param.replace(',','.')
			idx = int(uncertainty[uncertainty['Variables'] == param_aux].index[0])
			if type(uncertainty.loc[idx]['DEF_VALUES']) == str:
				inits = uncertainty.loc[idx]['DEF_VALUES'].split(',')
				mins = uncertainty.loc[idx]['MIN'].split(',')
				maxs = uncertainty.loc[idx]['MAX'].split(',')
				ranges = self.params_def([float(x) for x in mins], [float(x) for x in inits], [float(x) for x in maxs])
			else:
				inits = uncertainty.loc[idx]['DEF_VALUES']
				mins = uncertainty.loc[idx]['MIN']
				maxs = uncertainty.loc[idx]['MAX']
				ranges = self.params_def(mins, inits, maxs)
			ranges_to_plot = self.params_def(plot_params.loc[idx]['MIN'],plot_params.loc[idx]['DEF_VALUES'],plot_params.loc[idx]['MAX'])
			for i in range(len(self.experimental)):
				if float(self.experimental.loc[i][param]) >= -1  and float(self.experimental.loc[i][param]) <= 1:
					value = self.interpol(ranges,self.experimental.loc[i][param])
					full_param_value.append(str(value))
					value_to_plot = self.interpol(ranges_to_plot, self.experimental.loc[i][param])
					full_param_value_to_plot.append(float(value_to_plot[0]))
				else:
					full_param_value.append(str(self.experimental.loc[i][param]))
					full_param_value_to_plot.append(float(self.experimental.loc[i][param]))
			paramvalues[param] = full_param_value
			plot_paramvalues[param] = full_param_value_to_plot
		return paramvalues, plot_paramvalues

	def writeTXTFile(self, oFfolder, OF_values):
		self.sorted_OF_values = dict(sorted(OF_values.items(), key=lambda item: item[1][0]))
		with open(f'{self.simdir}\Sens_Results\{oFfolder}\{oFfolder}_sens_results.txt', 'w') as file:
			file.write(
				f'Minimum OF value obtained in {min(self.sorted_OF_values, key=self.sorted_OF_values.get)}. \nOF value = {min(self.sorted_OF_values.values())[0]}\n\n')
			for simu, result in self.sorted_OF_values.items():
				file.write(f'{simu}: {result[0]}\n')
				for well in self.OF.geowell.keys():
					wellprops = result[1].loc['Well OF Value',well]
					file.write(f'Well {well}:\t {wellprops}\n')

	def writePROOFstats(self, oFfolder, OF_values, evidence_combining_method):
		if os.path.isdir(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats') == False:
			os.makedirs(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats')

		for simu, result in OF_values.items():
			if os.path.isdir(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats\\Simu_{simu}') == False:
				os.makedirs(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats\\Simu_{simu}')
			print(f'Writing statistics for {simu}')
			stats = result[2]
			with pd.ExcelWriter(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats\\Simu_{simu}\\Well_Sequences.xlsx') as writer1:
				for well in stats:
					if evidence_combining_method =='belief':
						geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends,  df_prob_f], [T_Values, F_Values], full_well_belief, belief  = stats[well]
					elif evidence_combining_method =='mean':
						[geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends, df_prob_f], [T_Values, F_Values]] = stats[well]
					geowell.loc['ID', :] = 'GeoWell'
					simwell.loc['ID', :] = 'SimWell'
					welltable = pd.concat([geowell, simwell], axis=1)
					welltable.to_excel(writer1, sheet_name=f'{well}')

			with pd.ExcelWriter(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats\\Simu_{simu}\\Well_Statistics.xlsx') as writer2:
				for well in stats:
					if evidence_combining_method =='belief':
						geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends,  df_prob_f], [T_Values, F_Values], full_well_belief, belief  = stats[well]
					elif evidence_combining_method =='mean':
						[geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends, df_prob_f], [T_Values, F_Values]] = stats[well]
					geowellstats.loc['ID', :] = 'GeoWell'
					geowellstats.loc['-', :] = '-'
					simwellstats.loc['ID', :] = 'SimWell'
					simwellstats.loc['-', :] = '-'
					T_F_values = pd.concat([T_Values, F_Values], axis = 0)
					T_F_values.loc['-', :] = '-'
					df_prob_t.index = ['To_Be_Prob (Mean)']
					df_prob_f.index = ['To_Be_Prob_(Var)']
					df_prob_t_tends.index = ['To_Be_Prob_(Tends Mean)']
					probabilities = pd.concat([df_prob_t, df_prob_t_tends, df_prob_f], axis=0)
					probabilities.loc['-', :] = '-'
					probabilities.fillna('-', inplace=True)
					statstable = pd.concat([geowellstats, simwellstats, T_F_values, probabilities], axis=0)
					statstable.to_excel(writer2, sheet_name=f'{well}')

			if evidence_combining_method =='belief':
				with pd.ExcelWriter(f'{self.simdir}\\Sens_Results\\{oFfolder}\\stats\\Simu_{simu}\\Well_evidences_belief.xlsx') as writer3:
					for well in stats:
						geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends,  df_prob_f], [T_Values, F_Values], full_well_belief, belief = stats[well]
						full_well_belief.to_excel(writer3, sheet_name=f'{well}')

	def writeOFsTable(self,OF_values, oFfolder, param_table, param_table_to_plot):
		#TOTAL OF VALUES
		ofs_to_excel = pd.DataFrame({'simulations': [x for x, y in OF_values.items()], 'OF value': [y[0] for x, y in OF_values.items()]})
		#REGISTERING EXP DESIGN VALUES
		if self.expdesign != '':
			ofs_to_excel = pd.concat([ofs_to_excel, param_table_to_plot],axis=1)
			ofs_full_curves = pd.concat([ofs_to_excel, param_table],axis=1)
		else:
			ofs_full_curves = ofs_to_excel
		#SORTING VALUES TO MATCH
		ofs_full_curves = ofs_full_curves.sort_values(by="OF value")
		wellofs_results = {}
		with pd.ExcelWriter(f'{self.simdir}\Sens_Results\{oFfolder}\{oFfolder}_values_table_{self.case}.xlsx') as writer:
			ofs_full_curves.to_excel(writer, sheet_name= 'sorted_OF_Values_vs_Params')

			for well in self.OF.geowell.keys():
				wellofs = []
				for simu, result in OF_values.items():
						wellofs.append(result[1].loc['Well OF Value',well])
				wellofs_to_excel = pd.DataFrame({'simulations':  [x for x, y in OF_values.items()], f'Well{well}': wellofs})
				if self.expdesign != '':
					wellofs_to_excel = pd.concat([wellofs_to_excel, param_table_to_plot], axis=1)
				wellofs_to_excel = wellofs_to_excel.sort_values(by=f'Well{well}')
				wellofs_to_excel.to_excel(writer, sheet_name=f'Well{well}')

				wellofs_results[f'Well{well}'] = wellofs_to_excel
		return ofs_to_excel, wellofs_results

	def plotWellSimsOFs(self, oFfolder, param_table_to_plot, wellofs_results):
		for param in param_table_to_plot:
			fig = matplotlib.pyplot.gcf()
			fig.set_size_inches(12, 5)

			for well,result in wellofs_results.items():
				result = result.sort_values(by=f'{param}')
				plt.plot(result[param], result[f'{well}'], marker='o', linewidth=1.5, label = f'{well}')

			plt.legend(title='Well', title_fontsize=15, loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 'medium')
			plt.xlabel(f'{param}')
			plt.ylabel(f'Wells OF values')
			plt.title(f'Wells OF values versus {param}')
			plt.savefig(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}\WellOFs_vs{pathvalidate.sanitize_filename(param)}_values.png')
			plt.close()

	def processingSCOOFcontribution(self):
		#PROCESSING SCOOF CONTRIBUTIONS TO PLOT
		termscontribution = pd.DataFrame(np.zeros((len(self.OF.simulation.SCOOFoutputdata), len(self.OF.simulation.SCOOFheader))), columns=self.OF.simulation.SCOOFheader)
		for x in range(0,len(self.OF.simulation.SCOOFoutputdata)):
			termscontribution.iloc[x] = self.OF.simulation.SCOOFoutputdata[x]
		columns = [f'Well {x}' for x in self.OF.geowell.keys()]
		columns.append('Total')
		facies_terms_to_plot = pd.DataFrame(np.zeros((len(self.OF.simulation.SCOOFoutputdata), len(self.OF.geowell)+1)), columns = columns)
		thickness_terms_to_plot = pd.DataFrame(np.zeros((len(self.OF.simulation.SCOOFoutputdata), len(self.OF.geowell)+1)), columns = columns)
		for i in range(0,len(self.OF.simulation.SCOOFoutputdata)):
			for well in self.OF.geowell.keys():
				facies_terms_to_plot.loc[i].loc[f'Well {well}'] = termscontribution.loc[i][f'{well} - Facies term'] * termscontribution.loc[i]['Alpha']
				thickness_terms_to_plot.loc[i].loc[f'Well {well}'] = termscontribution.loc[i][f'{well} - Thickness term'] * termscontribution.loc[i][f'{well} - Beta']
			facies_terms_to_plot.loc[i].loc[f'Total'] = np.sum(facies_terms_to_plot.loc[i][:])/len(self.OF.geowell)
			thickness_terms_to_plot.loc[i].loc[f'Total'] = np.sum(thickness_terms_to_plot.loc[i][:])/len(self.OF.geowell)
		return facies_terms_to_plot, thickness_terms_to_plot

	def plotGlobalContributions(self, oFfolder, ofs_to_excel, param_table_to_plot, facies_terms_to_plot, thickness_terms_to_plot):
		facies_terms_to_plot = pd.concat([facies_terms_to_plot, param_table_to_plot, ofs_to_excel['simulations']],axis=1)
		for param in param_table_to_plot:
			facies_terms_to_plot = facies_terms_to_plot.sort_values(by=f'{param}')
			sims = []
			for sim in facies_terms_to_plot['simulations']:
				emp_str = ""
				for m in sim:
					if m.isdigit():
						emp_str = emp_str + m
				sims.append(int(emp_str))

			fig, ax = plt.subplots()
			fig.set_size_inches(12, 5)
			plt.plot(facies_terms_to_plot[param], facies_terms_to_plot['Total'], marker='o', linewidth = 1.5, markerfacecolor ='r', color='r')
			plt.xlabel(f'{param}')
			plt.ylabel('Facies Contribution')
			plt.title(f'Facies Global Contributions (Term * alpha) versus {param}')
			x = facies_terms_to_plot[param].reset_index(drop=True)
			y = facies_terms_to_plot['Total'].reset_index(drop=True)
			n = [str(t) for t in sims]
			for i, txt in enumerate(n):
				ax.annotate(txt, (x[i], y[i]),horizontalalignment='center', verticalalignment='top',xytext=(0,15), size = 10, textcoords='offset points')
			plt.savefig(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}\\{oFfolder}_Global_Facies_Contribution_vs_{pathvalidate.sanitize_filename(param)}.png')
			plt.close()

		thickness_terms_to_plot = pd.concat([thickness_terms_to_plot, param_table_to_plot, ofs_to_excel['simulations']], axis =1)
		for param in param_table_to_plot:
			thickness_terms_to_plot = thickness_terms_to_plot.sort_values(by=f'{param}')
			sims = []
			for sim in thickness_terms_to_plot['simulations']:
				emp_str = ""
				for m in sim:
					if m.isdigit():
						emp_str = emp_str + m
				sims.append(int(emp_str))

			fig, ax = plt.subplots()
			fig.set_size_inches(12, 5)
			plt.plot(thickness_terms_to_plot[param], thickness_terms_to_plot['Total'], marker='o', linewidth=1.5, markerfacecolor='b', color='b')
			plt.xlabel(f'{param}')
			plt.ylabel('Thickness Contribution')
			plt.title(f'Thickness Global Contributions (Term * beta) versus {param}')
			x = thickness_terms_to_plot[param].reset_index(drop=True)
			y = thickness_terms_to_plot['Total'].reset_index(drop=True)
			n = [str(t) for t in sims]
			for i, txt in enumerate(n):
				ax.annotate(txt, (x[i], y[i]),horizontalalignment='center', verticalalignment='top',xytext=(0,15), size = 10, textcoords='offset points')
			plt.savefig(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}\\{oFfolder}_Global_Thickness_Contribution_vs_{pathvalidate.sanitize_filename(param)}.png')
			plt.close()

	def plotWellsContributions(self, oFfolder, ofs_to_excel, param_table_to_plot, facies_terms_to_plot, thickness_terms_to_plot):
		facies_terms_to_plot = pd.concat([facies_terms_to_plot, param_table_to_plot, ofs_to_excel['simulations']], axis=1)
		facies_terms_to_plot.drop('Total',axis=1, inplace=True)
		facies_terms_to_plot.drop('simulations',axis=1,inplace=True)
		for param in param_table_to_plot:
			fig = matplotlib.pyplot.gcf()
			fig.set_size_inches(12, 5)
			facies_terms_to_plot = facies_terms_to_plot.sort_values(by=f'{param}')
			for well in facies_terms_to_plot.columns:
				if well != param:
					plt.plot(facies_terms_to_plot[param], facies_terms_to_plot[f'{well}'], marker='o', linewidth=1.5, label=f'{well}')

			plt.legend(title='Well', title_fontsize=15, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='medium')
			plt.xlabel(f'{param}')
			plt.ylabel(f'Wells Facies Contribution')
			plt.title(f'Wells Facies Contribution versus {param}')
			plt.savefig(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}\\{oFfolder}_Wells_Facies_Contribution_vs_{pathvalidate.sanitize_filename(param)}.png')
			plt.close()

		thickness_terms_to_plot = pd.concat([thickness_terms_to_plot, param_table_to_plot, ofs_to_excel['simulations']], axis=1)
		thickness_terms_to_plot.drop('Total',axis=1, inplace=True)
		thickness_terms_to_plot.drop('simulations',axis=1,inplace=True)
		for param in param_table_to_plot:
			fig = matplotlib.pyplot.gcf()
			fig.set_size_inches(12, 5)
			thickness_terms_to_plot = thickness_terms_to_plot.sort_values(by=f'{param}')
			for well in thickness_terms_to_plot.columns:
				if well != param:
					plt.plot(thickness_terms_to_plot[param], thickness_terms_to_plot[f'{well}'], marker='o', linewidth=1.5,label=f'{well}')

			plt.legend(title='Well', title_fontsize=15, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='medium')
			plt.xlabel(f'{param}')
			plt.ylabel(f'Wells Thickness Contribution')
			plt.title(f'Wells Thickness Contribution versus {param}')
			plt.savefig(f'{self.simdir}\\Sens_Results\\{oFfolder}\\{pathvalidate.sanitize_filename(param)}\\{oFfolder}_Wells_Thickness_Contribution_vs_{pathvalidate.sanitize_filename(param)}.png')
			plt.close()