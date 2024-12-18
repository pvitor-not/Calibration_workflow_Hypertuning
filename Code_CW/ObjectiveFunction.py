#Import libraries
from abc import ABC, abstractmethod
import Code_CW.Functions as Functions
import Code_CW.Optimizer as Optimizer
import Code_CW.Objects as Objects
import numpy as np
import math
import os
import pandas as pd
from scipy.stats.stats import _ttest_finish
from scipy.stats import ttest_ind
from scipy.stats import ttest_ind_from_stats
from scipy import stats

class Objective_function(ABC):
	def __init__(self):
		pass

	def set_objective_function(self, case, simdir: str, arcfilebase: str, origin: list, well_files: list,
							   extract_columns: list, facies_definition: object,simulation: object,
							   markers_table ='', litho_class_ids = '',carbo_ids = '', mud_ids='', sand_ids='', gravel_ids ='',
							   alpha='',beta='',internal_optimization_method='', evidence_combining_method = 'mean', total_OF_computing_method = 'mean') -> float:
		"""
		Parameters
		simdir: directory path of the input files.
		arcfilebase: .arc file used as base case.
		origin: origin coordinates.
		well_files: directory path of the well data .las files.
		extract_columns: columns names to be extracted from the DionisosFlow output files.
		facies_definition: object with information about the facies definition.
		simulation: object to save all the output data of the calibration (uncertain parameters and FO values).
		alpha:  weight for the facies distance (SCOOF only)
		beta:  weight for the thickness distance (SCOOF only)
		internal_optimization_method: internal optimization method (SCOOF only)

		Returns the total OF value.
		-------
		"""
		self.simdir = simdir
		self.arcfilebase = arcfilebase
		self.origin = origin
		self.well_files = well_files
		self.extract_columns = extract_columns
		self.facies_definition = facies_definition
		self.simulation = simulation
		self.set_alpha = alpha
		self.set_beta = beta
		self.allparametersprint = []
		self.markers_table = markers_table
		self.litho_class_ids = litho_class_ids
		self.carbo_ids = carbo_ids
		self.mud_ids = mud_ids
		self.sand_ids = sand_ids
		self.gravel_ids = gravel_ids
		self.case = case
		self.internal_optimization_method = internal_optimization_method
		self.evidence_combining_method = evidence_combining_method
		self.total_OF_method = total_OF_computing_method

	def set_optmization(self, uncertain_parameters, enginefolder: str, internal_optimizer: object,external_optimizer: object,
						reference_value = []):
		self.uncertain_parameters = uncertain_parameters
		self.reference_value = reference_value
		self.enginefolder = enginefolder
		self.internal_optimizer = internal_optimizer
		self.external_optimizer = external_optimizer

	def replacement(self, inversible_parameters, uncertain_parameters):
		"""
		Organize the values to be replaced in the uncertain parameters
		"""
		print('Inversible parameters are: ', inversible_parameters)
		outputpath = f'output{self.simulation.counter}'
		self.replaces = []
		for count, parameter in enumerate(self.uncertain_parameters.iloc[:, 0]):
			parameter_value = [parameter, inversible_parameters.item(count)]
			self.replaces.append(parameter_value)
		self.replaces.append(['OutputPath', f'{outputpath}'])
		self.arcinputfile, self.replaces = Functions.replaceInversibleParameters(self.replaces, self.simdir, self.arcfilebase, uncertain_parameters)
		return

	def new_simulation(self):
		"""
		Call the function to run Dionisos simulation
		"""
		Functions.runDionisosSim(self.simdir, self.enginefolder, self.arcinputfile)
		return

	def extract_simulation_data(self):
		# Object containing simulation data from base case.
		self.simdata = Functions.getSimData(self.simdir, self.arcinputfile, self.origin)
		return

	def extract_obs_well_data(self):
		self.geowell = Functions.getGeoWellData(self.well_files, Functions.getMeshCoordinates(self.simdata), self.extract_columns, self.markers_table)
		return

	def extract_sim_well_data(self):
		"""
		Call the functions to extract the results from Dionisos simulation
		"""
		# Dictionary contaning geowell objects getGeoWellData function.
		self.simwell_data = Functions.getSimWellData(self.well_files, self.simdir, Functions.getGeoWellRelativePosition(self.geowell),
									  Functions.toWellExtraction(self.simdir, self.arcinputfile, self.origin, self.extract_columns),
									  self.facies_definition)
		return

	def computeTotalOF_BOTH(self):
		"""
		Method to compute the total weighted OF value considering all geological wells

		Return a scalar (total OF value)
		"""
		wells_belief = pd.DataFrame(np.zeros((3,len(self.geowell.keys()))), columns = self.geowell.keys())
		wells_OF = pd.DataFrame(np.zeros((1, len(self.geowell.keys()))), columns=self.geowell.keys(), index = ['Well OF Value'])
		wellstats = {}
		for key in self.geowell.keys():
			self.simwell = self.simwell_data[key]  #Chamar informações do poço simulado no simulatedwell(key)
			self.well = self.geowell.get(key)

			if self.OF_type == 'S':
				self.beta = self.set_beta if self.set_beta != 0 else (1/self.well.totalthickness)
				self.simwell.distFacies(self.facies_definition)
				well_OF = self.computeWellOF()
				sequences = [self.well.sequence, self.simwell.table]
				wellstats[key] = sequences

			elif self.OF_type == 'P':
				if self.reference_value != []:
					self.well = self.ref_solution.get(key)
				well_OF, wellstats[key] = self.computeWellOF_fromsamples()
				if self.evidence_combining_method == 'belief':
					wells_belief[key] = wellstats[key][-1]
			wells_OF[key] = well_OF

		if self.OF_type =='S':
			total_OF_value = np.mean(wells_OF.iloc[0,:])
			print("Partial OF value: ", total_OF_value)
			return total_OF_value, wells_OF, wellstats

		elif self.OF_type == 'P':
			if self.total_OF_method =='belief':
				total_OF_value, total_OF_belief = self.BeliefFunctionTotalOFComputing(wells_belief)
				print("Partial OF value: ", total_OF_value)
				return total_OF_value, wells_OF, wellstats, total_OF_belief
			elif self.total_OF_method =='mean':
				total_OF_value = np.mean(wells_OF.loc['Well OF Value'])
				print("Partial OF value: ", total_OF_value)
				return total_OF_value, wells_OF, wellstats

	@abstractmethod
	def compute(self):
		"""
		Method to compute the total OF
		Returns the total OF value in float format
		"""
		pass

class Stratigraphic_OF(Objective_function):
	def __init__(self):
		self.OF_type = 'S'
		pass

	def computeSCOOF(self, cells):
		"""
		Compute the OF value for each candidate solution defined by the optimization method

		Parameters:
			cells: Array with candidate groupings for comparison with each simulated cell
		Return a scalar (OF value)
		"""

		if self.internal_optimizer == "AG" or self.internal_optimizer == "PSO":
			mi = 500
			cells = np.round(cells.astype(int), 0)
			if sum(cells) > last:
				return 1 + mi * abs(last - sum(cells))

		# cells = np.round(np.sort(cells).astype(int), 0)  # ascending order?? No!
		cells = np.round(cells.astype(int), 0)
		cells = tuple(cells)
		cells = (0,) + cells + (self.last,)  # combination
		OF = 0  # starts with FO value equal to zero

		for i in self.list_cells:  # counter with the number of simulated layers
			OF += (np.multiply(self.well.geocells_matrix[cells[i]], self.facies_simwell[i]).sum(axis=1) * self.alpha + abs(
				self.well.geothick[cells[i]] - self.simthick[i]) * self.beta).item(cells[i + 1] - cells[i])

		return OF

	def computeTermsContribution(self, cells):
		"""
		Compute the OF value for each term (facies and thickness)

		Return terms values
		"""
		cells = np.round(cells.astype(int), 0)
		cells = tuple(cells)
		OF = 0  # starts with FO value equal to zero
		facies_term = 0
		thickness_term = 0

		for i in self.list_cells:  # counter with the number of simulated layers
			# facies_term += (np.multiply(self.well.geocells_matrix[cells[i]], self.facies_simwell[i]).sum(axis=1) * self.alpha).item(cells[i + 1] - cells[i])
			# thickness_term += (abs(self.well.geothick[cells[i]] - self.simthick[i]) * self.beta).item(cells[i + 1] - cells[i])
			facies_term += (np.multiply(self.well.geocells_matrix[cells[i]], self.facies_simwell[i]).sum(axis=1)).item(cells[i + 1] - cells[i])
			if self.internal_optimization_method == 'DSR':
				# thickness_term += (abs(([j ** (1/2) for j in self.well.geothick[cells[i]]]) - (self.simthick[i]**(1/2)))).item(cells[i + 1] - cells[i])
				thickness_term += (abs(([j ** (1 / 2) for j in self.well.geothick[cells[i]]]) - (self.simthick[i] ** (1 / 2)))).item(cells[i + 1] - cells[i])
			else:
				thickness_term += (abs(self.well.geothick[cells[i]] - self.simthick[i])).item(cells[i + 1] - cells[i])

			# OF += (np.multiply(self.well.geocells_matrix[cells[i]], self.facies_simwell[i]).sum(axis=1) * self.alpha + abs(
			# 	([j ** (1/2) for j in self.well.geothick[cells[i]]]) - (self.simthick[i]**(1/2))) * self.beta).item(cells[i + 1] - cells[i])

		return facies_term, thickness_term

	def computeWellOF(self):
		"""
		Compute the OF value for each interpreted geological well

		Return a scalar (OF value)
		"""
		# interp_well = []
		self.facies_simwell = np.array(self.simwell.normdata[self.well.facies_list])
		self.simthick = np.array(self.simwell.normdata['Thickness'])
		self.list_cells = range(self.simwell.nb_cells)
		self.last = self.well.nb_data

		best_variable, min_OF = self.internal_optimizer.optimize(self.computeSCOOF, np.zeros(self.simwell.nb_cells - 1), np.ones(self.simwell.nb_cells - 1) * self.well.nb_data, 100, 100, self)

		facies_term, thickness_term = self.computeTermsContribution(best_variable)
		self.SCOOF_terms.append(facies_term)
		self.SCOOF_terms.append(thickness_term)
		self.SCOOF_terms.append(self.beta)

		cells = best_variable.astype(int)
		cells = np.sort(cells)
		cells = tuple(cells)
		cells = (0,) + cells + (self.last,)
		my_cells = cells
		well_thick = [self.well.geothick[my_cells[i]][my_cells[i + 1] - my_cells[i]] for i in range(self.simwell.nb_cells)]
		interpreted_well = [self.well.geocells_matrix[my_cells[i]][my_cells[i + 1] - my_cells[i]] for i in
							range(self.simwell.nb_cells)]
		interp_well = interpreted_well * np.array(well_thick)[:, None]
		interp_well = np.concatenate((interp_well, (interp_well.sum(axis=1)).reshape(-1, 1)), axis=1)
		print("OF value: ", min_OF)

		return min_OF

	def computeTotalOF(self):
		wells_OF = pd.DataFrame(np.zeros((1, len(self.geowell.keys()))), columns=self.geowell.keys(), index=['Well OF Value'])
		wellstats = {}
		for key in self.geowell.keys():
			self.simwell = self.simwell_data[key]
			self.well = self.geowell.get(key)
			self.beta = self.set_beta if self.set_beta != 0 else (1 / self.well.totalthickness)

			self.simwell.distFacies(self.facies_definition)

			wells_OF[key] = self.computeWellOF()
			wellstats[key]  = [self.well.sequence, self.simwell.table]

		total_OF_value = np.mean(wells_OF.iloc[0, :])
		print("Partial OF value: ", total_OF_value)
		return total_OF_value, wells_OF, wellstats

	def compute(self, inversible_parameters):
		"""
		Method to compute the total OF

		Returns the total OF value in float format
		"""
		self.SCOOF_terms = []
		self.allparametersprint.append(inversible_parameters)
		self.replacement(inversible_parameters,self.uncertain_parameters) #create new DF simulation based on realization values
		self.new_simulation() #run DF simulation
		self.extract_simulation_data()
		self.extract_obs_well_data()
		self.alpha = self.set_alpha if self.set_alpha != 0 else (1 / np.sqrt(self.facies_definition.nb_prop))
		self.extract_sim_well_data()  #extract DF results
		total_OF_value, wells_OF, wellstats = self.computeTotalOF() # compute OF

		self.SCOOF_terms.append(self.alpha)
		self.simulation.saveSCOOFterms(self.SCOOF_terms)
		self.simulation.saveoutput(self.replaces, total_OF_value)#save results of each round
		self.simulation.writeSCOOFtoafile(self.simdir, self.OF_type, self.geowell)
		self.simulation.writeTotalOFresults(self.simdir, self.OF_type, total_OF_value, wells_OF)
		self.simulation.writeSCOOFfacies(self.simdir, self.OF_type, wellstats)
		self.simulation.writetoafile(self.simdir,self.OF_type)
		self.simulation.updatecounter()
		print('Total OF value is: ', total_OF_value)

		print(self.allparametersprint)

		return (total_OF_value)

class Probabilistic_OF(Objective_function):
	def __init__(self):
		self.OF_type = 'P'
		pass

	def set_lithologies(self):
		self.litho_ids = list(self.litho_class_ids.keys())
		self.litho_classes = list(self.litho_class_ids.values())
		self.carbo_classes = [self.litho_class_ids[f'{x}'] for x in self.carbo_ids]
		self.terr_ids = [x for x in self.litho_ids if x not in self.carbo_ids]
		self.terr_classes = [self.litho_class_ids[f'{x}'] for x in self.terr_ids]
		self.mud_classes = [self.litho_class_ids[f'{x}'] for x in self.mud_ids]
		self.sand_classes = [self.litho_class_ids[f'{x}'] for x in self.sand_ids]
		self.gravel_classes = [self.litho_class_ids[f'{x}'] for x in self.gravel_ids]
		pass

	def ftest(self,f, N):
		'''
		Function to execute f test
        Parameters:
            v1: variance 1
            v2: variance 2
            n1: sample size 1
            n2: sample size 2
        Returns
            p: p-value of F test statistic
        '''
		dfn = N - 2  # define degrees of freedom numerator
		dfd = N - 2  # define degrees of freedom denominator
		p = 2*(1 - stats.f.cdf(f, dfn, dfd))  #find p-value of F test statistic
		return p

	def applyThicknessSqRoot(self, well):
		well['Thickness'] = np.sqrt(well['Thickness'])
		return well

	def applyProbabilitiesLimits(self, probdata):
		probabilities = probdata.copy(deep=True)
		for x in probabilities.columns:
			if float(probabilities[x]) > 0.9:
				probabilities[x] = 0.9
			elif float(probabilities[x]) < 0.1:
				probabilities[x] = 0.1
		return probabilities

	def computeT_F_values(self, geostats, simstats):
		T_values = pd.DataFrame(np.zeros((1, len(geostats.columns))), columns=geostats.columns, index = ['T_Value'])
		F_values = pd.DataFrame(np.zeros((1, len(geostats.columns))), columns=geostats.columns, index = ['F_Value'])
		for property in geostats.columns:
			#Calculating T_Values
			T_values.loc['T_Value', property] = abs(geostats.loc['Mean', property] - simstats.loc['Mean', property]) / np.sqrt(
				(np.mean([geostats.loc['Var', property], simstats.loc['Var', property]])/ (self.N - 1)))

			if math.isnan(T_values.loc['T_Value', property]):
				if (geostats.loc['Mean', property] == 0 and simstats.loc['Mean', property] == 0):
					T_values.loc['T_Value', property] = 0
				else:
					T_values.loc['T_Value', property] = 999
			#Calculating F_Values
			condition = simstats.loc['Var', property] / geostats.loc['Var', property]
			if condition.item() >=1 and condition.item() != math.inf:
				F_values.loc['F_Value', property] = condition
			elif condition.item() <1 and condition.item() >0:
				F_values.loc['F_Value', property] = geostats.loc['Var', property] / simstats.loc['Var', property]
			else:
				F_values.loc['F_Value', property] = 0
				print("divisão por zero na condition")

		return T_values, F_values

	def concatenateprobabilities(self, t_probs, f_probs, trend_probs):
		t = pd.DataFrame(t_probs.values, columns = [f'{x} (Mean)' for x in t_probs.columns], index = ['To_Be_Prob'])
		f = pd.DataFrame(f_probs.values, columns=[f'{x} (Var)' for x in f_probs.columns], index=['To_Be_Prob'])
		trend = pd.DataFrame(trend_probs.values, columns=[f'{x} (Trend Mean)' for x in trend_probs.columns], index=['To_Be_Prob'])
		df_dist = pd.concat([t, f, trend], axis=1)
		return df_dist

	def ttestfromsamples(self, geowell, simwell, interest_properties):
		prob_values = pd.DataFrame(np.zeros((1, len(interest_properties))), columns = interest_properties, index= ['To_Be_Prob'])
		for property in interest_properties:
			prob_values.loc['To_Be_Prob', property] = ttest_ind(geowell[property],simwell[property],equal_var=False,alternative= 'two-sided')[1]
		return prob_values

	def ttestfromstats(self,T_inputs, interest_properties):
		prob_values = pd.DataFrame(np.zeros((1, len(interest_properties))), columns = interest_properties, index= ['To_Be_Prob'])
		for property in interest_properties:
			prob_values.loc['To_Be_Prob', property] = _ttest_finish(self.N-2,T_inputs.loc['T_Value', property],alternative='two-sided')[1]
		return prob_values

	def ftestfromstats(self,F_inputs, interest_properties):
		prob_values = pd.DataFrame(np.zeros((1, len(interest_properties))), columns= interest_properties, index=['To_Be_Prob'])

		for property in interest_properties:
			prob_values.loc['To_Be_Prob', property] = self.ftest(F_inputs.loc['F_Value', property],self.N)
		return prob_values

	def set_reference_solution(self):
		ref_replace = self.reference_value
		ref_replace.append(['OutputPath', 'outputR'])
		ref_arcinput = Functions.replaceInversibleParameters(ref_replace, self.simdir, self.arcfilebase, self.uncertain_parameters)[0]
		Functions.runDionisosSim(self.simdir, self.enginefolder, ref_arcinput)
		self.ref_solution = Functions.getSimWellData(self.well_files, self.simdir, Functions.getGeoWellRelativePosition(self.geowell),
									  Functions.toWellExtraction(self.simdir, ref_arcinput, self.origin, self.extract_columns),
									  self.facies_definition)
		return

	def getDominantLitho(self,df_input):
		dom_litho = []
		columns = self.litho_classes
		for i in range(len(df_input)):
			layer_litho = df_input.loc[i, columns]
			if sum(layer_litho) != 0:
				layer_litho = layer_litho.idxmax()
				dom_litho.append(layer_litho)
			else:
				dom_litho.append('0')
		for id, litho in self.litho_class_ids.items():
			dom_litho = list(map(lambda x: x.replace(litho, id), dom_litho))
		dom_litho = [int(x) for x in dom_litho]
		df_input = df_input.drop(self.litho_classes,axis = 1)
		df_input['Litho_Dom'] = dom_litho

		return df_input

	def setTextClass(self,df_input):
		for i in range(0,len(df_input)):
			mud_idx = 1
			mud = np.sum(df_input.loc[i,self.mud_classes])
			sand_idx=3
			sand = np.sum(df_input.loc[i,self.sand_classes])
			gravel_idx=5
			gravel = np.sum(df_input.loc[i].loc[self.gravel_classes])

			if (mud>sand and mud>gravel) or (mud >= 0.999):
				id = mud_idx
			elif(sand>=mud and sand>gravel) or (sand>=0.999):
				id = sand_idx
			elif(gravel>mud and gravel>sand) or (gravel>=0.999):
				id = gravel_idx
			elif (mud ==0 and sand ==0 and gravel ==0):
				id = 0
			df_input.loc[i, 'Ltx'] = id

			if id <2:
				if df_input.loc[i,'Terrigenos'] >0.5:
					df_input.loc[i,'Litho'] = 1
				else:
						df_input.loc[i, 'Litho'] = 2
			elif (id>=2 and id<4):
				if df_input.loc[i,'Terrigenos'] >0.5:
					df_input.loc[i,'Litho'] = 3
				else:
					df_input.loc[i, 'Litho'] = 4
			elif id >=4:
				if df_input.loc[i,'Terrigenos'] >0.5:
					df_input.loc[i,'Litho'] = 5
				else:
					df_input.loc[i, 'Litho'] = 6
		return df_input

	def SimWellPreProcess(self,well):
		"""
            Prepare the data from the simulator

            Parameters:
                df_input: Simulated_Well Object
            Return the statistics results from the simulated properties
        """
		interest_properties = list(self.litho_classes)
		interest_properties.append('Bathymetry')
		new_df_sim = pd.DataFrame(data=None)
		dfaux = well.rawdata

		for property in dfaux.columns:
				if property in interest_properties:
					new_df_sim[property] = dfaux[property]
		new_df_sim['Thickness'] = dfaux['StratigraphicThickness']
		new_df_sim.reset_index(inplace=True, drop=True)

		new_df_sim.assign(Bat_Faixa=0)  # Adding a new column to the dataframe
		# Calculate the bathymetry using the definition
		bathymetryRange = [50, 20, 5 ,-1, -200]
		faixa = [4, 3, 2, 1, 0]

		#Setting bathymetry range for bathymetry value of each layer
		condition = new_df_sim.Bathymetry.reset_index(drop=True)
		for i in range(len(condition)):
			for j in range(len(bathymetryRange)):
				if condition.loc[i] > bathymetryRange[j]:
					new_df_sim.loc[i, 'Bat_Faixa'] = faixa[j]
					break
				elif condition.loc[i] < bathymetryRange[len(bathymetryRange) - 1]:
					new_df_sim.loc[i, 'Bat_Faixa'] = faixa[len(faixa) - 1]
		# Calculate the Terrigenos and Clastos Grossos
		new_df_sim.assign(Terrigenos=0)
		for i in range(len(new_df_sim.Bathymetry)):
			if (np.sum(new_df_sim.loc[i, self.litho_classes])) != 0:
				new_df_sim.loc[i, 'Terrigenos'] = np.sum(new_df_sim.loc[i, self.terr_classes]) / (
					np.sum(new_df_sim.loc[i, self.litho_classes]))
			else:
				new_df_sim.loc[i, 'Terrigenos'] = 0.0

		return (new_df_sim)

	def GeoWellPreProcess(self, well):
		"""
                Prepare the data from obs well

                Parameters:
                    df_input: Obs_Well Object
                    N: the size of the simulated data
                Return the statistics results from the obs properties
        """
		new_df_poço = pd.DataFrame.reset_index(well.sequence)
		# Compute the accumulated thickness
		for i in range(len(well.sequence) - 1, -1, -1):
			if i == len(well.sequence) - 1:
				new_df_poço.loc[i, 'esp_a'] = 0
			elif i == len(well.sequence) - 2:
				new_df_poço.loc[i, 'esp_a'] = new_df_poço.loc[i + 1, 'esp_a'] + new_df_poço.loc[i, 'Thickness'] + \
											  new_df_poço.loc[i + 1, 'Thickness']
			else:
				new_df_poço.loc[i, 'esp_a'] = new_df_poço.loc[i + 1, 'esp_a'] + new_df_poço.loc[i, 'Thickness']

		# Dividing the DataFrame in the same number of layers from the simulated data
		div = self.N + 1
		if len(new_df_poço.lithology) != div:
			columns = ['Esp_Max', 'Esp_Min', 'Intervalo']
			columns.extend(self.litho_ids)
			df2 = pd.DataFrame(0, index=np.array(range(div)), columns= columns)

			k = new_df_poço.loc[0, 'esp_a'] / div
			for i in range(div):
				df2.loc[i, 'Esp_Max'] = k * (div - i)
				df2.loc[i, 'Esp_Min'] = k * (div - i - 1)
				df2.loc[i, 'Intervalo'] = div - i

		# Compute the proportion of litofacies per thickness layer
		c = 0
		j = 0
		soma = 0
		while c < len(new_df_poço.lithology) - 1:
			c += 1
			if new_df_poço.loc[c, 'lithology'] != new_df_poço.loc[c - 1, 'lithology']:
				v = new_df_poço.loc[c - 1, 'esp_a']
				while j < len(df2.Intervalo):
					if df2.Esp_Min[j] <= v <= df2.Esp_Max[j]:
						if df2.loc[j, str(int(new_df_poço.loc[c - 1, 'lithology']))] != 0:
							df2.loc[j, str(int(new_df_poço.loc[c - 1, 'lithology']))] = 0

						soma = sum(df2.loc[j, self.litho_ids])
						df2.loc[j, str(int(new_df_poço.loc[c - 1, 'lithology']))] += df2.Esp_Max[j] - v - soma
						df2.loc[j, str(int(new_df_poço.loc[c, 'lithology']))] += v - df2.Esp_Min[j]
						j += 1
						break
					elif v > df2.Esp_Max[j]:
						j -= 1
					else:
						df2.loc[j, str(int(new_df_poço.loc[c - 1, 'lithology']))] = k
						j += 1

		# Filling in the lasts thickness layers
		while j < len(df2.Intervalo):
			df2.loc[j, str(int(new_df_poço.loc[c, 'lithology']))] = k
			j += 1
		df3 = pd.DataFrame(0, index=np.array(range(div)),
						   columns=['EspMean', 'LitoMean', 'VarLito', 'BatMean', 'VarBat', 'TerrMean', 'VarTerr'])

		# Compute the mean of Bathymetry per thickness layer
		c = 0
		j = 0
		while c < len(new_df_poço.lithology) - 1:
			c += 1
			if new_df_poço.loc[c, 'bathymetry'] != new_df_poço.loc[c - 1, 'bathymetry']:
				v = new_df_poço.loc[c - 1, 'esp_a']
				while j < len(df2.Intervalo):
					if df2.Esp_Min[j] <= v <= df2.Esp_Max[j]:
						if df2.loc[j, str(int(new_df_poço.loc[c - 1, 'lithology']))] > df2.loc[
							j, str(int(new_df_poço.loc[c, 'lithology']))]:  # calcula batimetria média para o intervalo
							df3.loc[j, 'BatMean'] = new_df_poço.loc[c - 1, 'bathymetry']
						else:
							df3.loc[j, 'BatMean'] = new_df_poço.loc[c, 'bathymetry']
						j += 1
						break
					elif v > df2.Esp_Max[j]:
						j -= 1
					else:
						df3.loc[j, 'BatMean'] = new_df_poço.loc[c - 1, 'bathymetry']
						j += 1 % 4  #
		# Filling in the lasts thickness layers
		while j < len(df2.Intervalo):
			df3.loc[j, 'BatMean'] = new_df_poço.loc[c, 'bathymetry']
			j += 1
		# Compute the properties mean

		# Create a new DataFrame (Until in test for the obs_well data)
		df4 = pd.DataFrame(df2.loc[:, self.litho_ids])
		for i in range(div):
			if np.sum(df4.loc[i, self.litho_ids]) > 0:
				df4.loc[i, 'Terrigenos'] = np.sum(df4.loc[i, self.terr_ids])/(np.sum(df4.loc[i, self.litho_ids]))
			else:
				df4.loc[i, 'Terrigenos'] = 0

		for id, litho in self.litho_class_ids.items():
			df4.rename(columns={f'{id}': f'{litho}'}, inplace=True)
			df4[f'{litho}'] = df4[f'{litho}']/k

		df4['Bat_Faixa'] = df3.loc[:, 'BatMean']
		df4['Thickness'] = k
		return (df4)

	def StatisticalCalculation(self, df_input):
		df_stats = pd.DataFrame(data=[])

		#Set Textural Classes and Lithologies based on sediment proportions
		df_input = self.setTextClass(df_input)
		#Calculate Litho Stats (Based on Textural Classification)
		df_stats.assign(Ltx=0)
		df_stats.loc['Mean', 'Ltx'] = np.sum(df_input.loc[:,'Ltx']*df_input.loc[:,'Thickness'])/np.sum(df_input.loc[:,'Thickness'])

		int_var = 1
		min_lit_var = 0.8 * df_stats.loc['Mean', 'Ltx']
		lit_var = int_var + (np.var(df_input.loc[:,'Ltx'] * df_input.loc[:, 'Thickness'] / np.sum(df_input.loc[:, 'Thickness']),ddof=1))
		if lit_var < min_lit_var:
			df_stats.loc['Var', 'Ltx'] = min_lit_var
		else:
			df_stats.loc['Var', 'Ltx'] = lit_var

		#Calculate Terrigen Stats
		df_stats.assign(Terrigenos=0)
		df_stats.loc['Mean', 'Terrigenos'] = np.sum(df_input.loc[:, 'Terrigenos'] * df_input.loc[:, 'Thickness']) / np.sum(df_input.loc[:, 'Thickness'])

		min_terr_var = 0.1
		terr_var = df_stats.loc['Mean', 'Terrigenos'] *(1-df_stats.loc['Mean', 'Terrigenos'])
		if terr_var < min_terr_var:
			df_stats.loc['Var', 'Terrigenos'] = min_terr_var
		else:
			df_stats.loc['Var', 'Terrigenos'] = terr_var

		#Calculate Bathymetry Stats
		df_stats.assign(Bat_Faixa=0)
		df_stats.loc['Mean', 'Bat_Faixa'] = np.mean(df_input.loc[:, 'Bat_Faixa'])

		min_bat_var = 0.8 * df_stats.loc['Mean', 'Bat_Faixa']
		bat_var = np.var(df_input.loc[:, 'Bat_Faixa'],ddof=1)
		if bat_var < min_bat_var:
			df_stats.loc['Var', 'Bat_Faixa'] = min_bat_var
		else:
			df_stats.loc['Var', 'Bat_Faixa'] = bat_var

		#Calculate Thicknes Stats (Based on Square Root of Thicknesses)
		df_stats.assign(Thickness=0)
		df_stats.loc['Mean', 'Thickness'] = np.mean(np.sqrt(df_input.loc[:, 'Thickness']))

		min_thick_var = 0.2
		thick_var = np.var(np.sqrt(df_input.loc[:, 'Thickness']), ddof = 1)
		if thick_var <= min_thick_var:
			df_stats.loc['Var', 'Thickness'] = min_thick_var
		else:
			df_stats.loc['Var', 'Thickness'] = thick_var

		#Calculate Litho/Bat Variability
		nt = 0
		for i in range(0,len(df_input)-1):
			if df_input.loc[i,'Litho'] != df_input.loc[i+1,'Litho']:
				nt += 1
			else:
				if df_input.loc[i,'Bat_Faixa'] != df_input.loc[i+1, 'Bat_Faixa']:
					nt += 1

		df_stats.assign(TransFreq=0)
		df_stats.loc['Mean','TransFreq'] = nt/(len(df_input)-1)
		df_stats.loc['Var','TransFreq'] = df_stats.loc['Mean','TransFreq']/(1-df_stats.loc['Mean','TransFreq'])
		#Calculate Trend Stats (Based on Textural Classification divided by superior and inferior zones)
		n = int(len(df_input.Thickness)/2)
		intervalo = np.linspace(1, 0, n, dtype=float)

		df_stats.assign(Sup_Ltx_Trend=0)
		suplit = df_input.loc[0:n-1,'Ltx']
		p_tend_sup = stats.linregress(intervalo,suplit.iloc[::-1])[0:2]
		slope = p_tend_sup[0]
		intercept = p_tend_sup[1]
		df_stats.loc['Mean', 'Sup_Ltx_Trend'] = slope
		## Parameter 2 to storage the linregress values (slope and interception)
		df_stats.loc['Var', 'Sup_Ltx_Trend'] = 0
		for i in range(n - 1):
			df_stats.loc['Var', 'Sup_Ltx_Trend'] += ((intercept + slope * intervalo[i]) - df_input.loc[i,'Ltx']) ** 2 \
													/ len(intervalo)

		df_stats.assign(Inf_Ltx_Trend=0)
		inflit = df_input.loc[n:len(df_input.Thickness), 'Ltx']
		p_tend_inf = stats.linregress(intervalo, inflit.iloc[::-1])[0:2]
		slope = p_tend_inf[0]
		intercept = p_tend_inf[1]
		df_stats.loc['Mean', 'Inf_Ltx_Trend'] = slope
		## Parameter 2 to storage the linregress values (slope and interception)
		df_stats.loc['Var', 'Inf_Ltx_Trend'] = 0
		for i in range(n - 1):
			df_stats.loc['Var', 'Inf_Ltx_Trend'] += ((intercept + slope * intervalo[i]) - df_input.loc[i, 'Ltx']) ** 2 \
														/ len(intervalo)

		return df_input, df_stats

	def PreProcessing(self):
		self.N = len(self.simwell.rawdata) - 1
		processed_simwell = self.SimWellPreProcess(self.simwell)

		if self.reference_value != []:
			processed_geowell = self.SimWellPreProcess(self.well)
		else:
			processed_geowell = self.GeoWellPreProcess(self.well)

		simwell, simwellstats = self.StatisticalCalculation(processed_simwell)
		geowell, geowellstats = self.StatisticalCalculation(processed_geowell)
		simwell_sqroot_samples = self.applyThicknessSqRoot(simwell)
		geowell_sqroot_samples = self.applyThicknessSqRoot(geowell)

		if self.case =='reference_solution_FO_zero':
			return simwell, simwellstats, simwell_sqroot_samples, simwell, simwellstats, simwell_sqroot_samples
		else:
			return geowell, geowellstats, geowell_sqroot_samples, simwell, simwellstats, simwell_sqroot_samples

	def BeliefFunctionEvidenceComputing(self, probabilities):
		belief_df = pd.DataFrame(np.zeros((3, len(probabilities.columns))), columns=probabilities.columns)
		belief_df = pd.concat([probabilities,belief_df], axis=0)
		belief_df.index = ['P(U|Ei)', 'P(U|E)', 'P(nU|E)', 'P(Ua|E)']
		belief_df.loc['P(U|E)'].iloc[0] = belief_df.loc['P(U|Ei)'].iloc[0]
		belief_df.loc['P(nU|E)'].iloc[0] = 1 - belief_df.loc['P(U|E)'].iloc[0]
		belief_df.loc['P(Ua|E)'].iloc[0] = 1 - (belief_df.loc['P(U|E)'].iloc[0]+belief_df.loc['P(nU|E)'].iloc[0])

		for i in range(1,len(belief_df.columns)):
			belief_df.loc['P(U|E)'].iloc[i] = belief_df.loc['P(U|Ei)'].iloc[i]*(belief_df.loc['P(U|E)'].iloc[i-1]+belief_df.loc['P(Ua|E)'].iloc[i-1])
			belief_df.loc['P(nU|E)'].iloc[i] = (1-belief_df.loc['P(U|Ei)'].iloc[i])*(belief_df.loc['P(nU|E)'].iloc[i-1]+belief_df.loc['P(Ua|E)'].iloc[i-1])
			belief_df.loc['P(Ua|E)'].iloc[i] = (belief_df.loc['P(U|Ei)'].iloc[i]*belief_df.loc['P(nU|E)'].iloc[i-1])+((1-belief_df.loc['P(U|Ei)'].iloc[i])*belief_df.loc['P(U|E)'].iloc[i-1])
		pUE = belief_df.iloc[1][-1]
		pnUE = belief_df.iloc[2][-1]
		pUaE = belief_df.iloc[3][-1]
		OF_value = pnUE
		probabilities = [pUE, pnUE, pUaE]
		full_well_belief = belief_df
		return OF_value, probabilities, full_well_belief

	def BeliefFunctionTotalOFComputing(self, wellsbelief):
		combination_columns = []
		for i in range(len(wellsbelief.columns)):
			if i <2:
				combination_columns.append(f'Well{wellsbelief.columns[i]}')
			else:
				combination_columns.append(f'c{i-1}')
				combination_columns.append(f'Well{wellsbelief.columns[i]}')
		combination_columns.append(f'c{len(wellsbelief.columns)-1}')

		fullbelief_df = pd.DataFrame(np.zeros((3, len(combination_columns))), columns=combination_columns)
		for well in wellsbelief:
			fullbelief_df[f'Well{well}'] = wellsbelief[well]

		fullbelief_df.index = ['P(U|E)', 'P(nU|E)', 'P(Ua|E)']
		for i in range(1,len(wellsbelief.columns)):
			t = fullbelief_df.columns.get_loc(f'c{i}')
			fullbelief_df.loc['P(U|E)', f'c{i}'] = (fullbelief_df.loc['P(U|E)'].iloc[t-1] * fullbelief_df.loc['P(U|E)'].iloc[t-2]) + (fullbelief_df.loc['P(U|E)'].iloc[t-1] * fullbelief_df.loc['P(Ua|E)'].iloc[t-2]) + (fullbelief_df.loc['P(U|E)'].iloc[t-2] * fullbelief_df.loc['P(Ua|E)'].iloc[t-1])
			fullbelief_df.loc['P(nU|E)', f'c{i}'] = (fullbelief_df.loc['P(nU|E)'].iloc[t-1] * fullbelief_df.loc['P(nU|E)'].iloc[t-2]) + (fullbelief_df.loc['P(nU|E)'].iloc[t-2] * fullbelief_df.loc['P(Ua|E)'].iloc[t-1]) + (fullbelief_df.loc['P(nU|E)'].iloc[t-1] * fullbelief_df.loc['P(Ua|E)'].iloc[t-2])
			fullbelief_df.loc['P(Ua|E)', f'c{i}'] = (fullbelief_df.loc['P(U|E)'].iloc[t-1] * fullbelief_df.loc['P(nU|E)'].iloc[t-2]) + (fullbelief_df.loc['P(U|E)'].iloc[t-2] * fullbelief_df.loc['P(nU|E)'].iloc[t-1]) + (fullbelief_df.loc['P(Ua|E)'].iloc[t-1]* fullbelief_df.loc['P(Ua|E)'].iloc[t-2])
		totalOF = fullbelief_df.loc['P(nU|E)'].iloc[-1]
		return totalOF, fullbelief_df

	def combineEvidences(self, method, df_probs=[], kwargs = []):
		geowell, geowellstats, simwell, simwellstats, T_Values, F_Values = kwargs
		df_prob_t, df_prob_f, df_prob_t_tends = df_probs
		df_dist = self.concatenateprobabilities(df_prob_t, df_prob_f, df_prob_t_tends)
		df_probabilities = self.applyProbabilitiesLimits(df_dist)
		df_probabilities.drop('Thickness (Var)', axis=1, inplace=True)

		if method =='belief':
			OF_value, belief, full_well_belief = self.BeliefFunctionEvidenceComputing(df_probabilities)
			wellstats = [geowell, geowellstats, simwell, simwellstats, [df_prob_t, df_prob_t_tends, df_prob_f],[T_Values, F_Values], full_well_belief, belief]

		elif self.evidence_combining_method =='mean':
			prob = np.mean(df_probabilities, axis=1)
			OF_value = 1-prob
			df_prob_t_lmt = self.applyProbabilitiesLimits(df_prob_t)
			df_prob_t_tends_lmt = self.applyProbabilitiesLimits(df_prob_t_tends)
			df_prob_f_lmt = self.applyProbabilitiesLimits(df_prob_f)
			wellstats = [geowell, geowellstats, simwell, simwellstats, [df_prob_t_lmt, df_prob_t_tends_lmt, df_prob_f_lmt],[T_Values, F_Values]]

		return OF_value, wellstats

	def computeWellOF_fromsamples(self):
		geowell, geowellstats, geowell_sqroot_samples, simwell, simwellstats, simwell_sqroot_samples = self.PreProcessing()

		T_Values, F_Values = self.computeT_F_values(geowellstats, simwellstats)
		df_prob_t = self.ttestfromsamples(geowell_sqroot_samples,simwell_sqroot_samples, ['Ltx', 'Bat_Faixa', 'Terrigenos', 'Thickness'])
		df_prob_t_tends = self.ttestfromstats(T_Values, ['Sup_Ltx_Trend', 'Inf_Ltx_Trend','TransFreq'])
		#df_prob_t_tends = self.ttestfromstats(T_Values, ['Sup_Ltx_Trend', 'Inf_Ltx_Trend'])
		df_prob_f = self.ftestfromstats(F_Values, ['Ltx', 'Bat_Faixa', 'Terrigenos', 'Thickness'])

		OF_value, wellstats = self.combineEvidences(self.evidence_combining_method, [df_prob_t, df_prob_f, df_prob_t_tends],
													kwargs = [geowell,geowellstats, simwell, simwellstats, T_Values, F_Values])

		OF_value = float(OF_value)
		print("OF value: ", OF_value)
		return OF_value, wellstats

	def computeTotalOF(self):
		wells_belief = pd.DataFrame(np.zeros((3, len(self.geowell.keys()))), columns=self.geowell.keys())
		wells_OF = pd.DataFrame(np.zeros((1, len(self.geowell.keys()))), columns=self.geowell.keys(), index=['Well OF Value'])
		wellstats = {}
		for key in self.geowell.keys():
			self.simwell = self.simwell_data[key]  # Chamar informações do poço simulado no simulatedwell(key)
			self.well = self.geowell.get(key)
			if self.reference_value != []:
				self.well = self.ref_solution.get(key)

			wells_OF[key], wellstats[key] = self.computeWellOF_fromsamples()
			if self.evidence_combining_method == 'belief':
				wells_belief[key] = wellstats[key][-1]
		if (self.evidence_combining_method =='mean' and self.total_OF_method=='belief'):
			wells_OF_to_belief = wells_OF.copy(deep=True)
			wells_OF_to_belief.loc['Well OF Value', :] = 1 - wells_OF_to_belief.loc['Well OF Value', :]
			total_OF_value, probabilities, full_well_belief = self.BeliefFunctionEvidenceComputing(wells_OF_to_belief)
			print("Partial OF value: ", total_OF_value)
			return total_OF_value, wells_OF, wellstats, full_well_belief

		elif (self.evidence_combining_method =='mean' and self.total_OF_method=='mean'):
			total_OF_value = np.mean(wells_OF.loc['Well OF Value'])
			print("Partial OF value: ", total_OF_value)
			return total_OF_value, wells_OF, wellstats

		elif (self.evidence_combining_method =='belief' and self.total_OF_method=='belief'):
			total_OF_value, total_OF_belief = self.BeliefFunctionTotalOFComputing(wells_belief)
			print("Partial OF value: ", total_OF_value)
			return total_OF_value, wells_OF, wellstats, total_OF_belief

	def compute(self, inversible_parameters):
		"""
		Method to compute the total OF

		Returns the total OF value in float format
		"""
		self.set_lithologies()
		self.allparametersprint.append(inversible_parameters)
		self.replacement(inversible_parameters, self.uncertain_parameters)
		self.new_simulation()
		self.extract_simulation_data()
		self.extract_obs_well_data()
		self.extract_sim_well_data()

		if self.reference_value != [] and self.simulation.counter == 0:
			print(f'!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!\nYou must provide uncertain parameters reference values to test MixedModel calibration using Probabilistic Objective Function (PROOF) !!!WARNING!!!')
			self.set_reference_solution()
			print(f'Reference Solution has been succesfully set.')
		self.lithologies = self.simdata.sediment_classes()

		if self.total_OF_method == 'belief':
			total_OF_value, wells_OF, fullstats, total_OF_belief = self.computeTotalOF()
			self.simulation.writeTotalBeliefresults(self.simdir, self.OF_type, total_OF_belief)
		elif self.total_OF_method =='mean':
			total_OF_value, wells_OF, fullstats = self.computeTotalOF()
			total_OF_belief = ''

		self.simulation.saveoutput(self.replaces, total_OF_value)
		self.simulation.writeTotalOFresults(self.simdir,self.OF_type,total_OF_value, wells_OF)
		self.simulation.writePROOFStats(self.simdir, self.OF_type, fullstats, self.evidence_combining_method)
		self.simulation.writetoafile(self.simdir,self.OF_type,kwargs=total_OF_belief)
		self.simulation.updatecounter()
		print('Total OF value is: ', total_OF_value)
		return (total_OF_value)