#Import libraries
from abc import ABC, abstractmethod
import Code_CW.Functions as Functions
import Code_CW.ObjectiveFunction as ObjectiveFunction
from pyswarm import pso
import numpy as np
from geneticalgorithm import geneticalgorithm as ga
from scipy.sparse import coo_matrix
from Code_CW import Dijkstra
from scipy.sparse.csgraph import dijkstra
import itertools

class Optimizer(ABC):
	def __init__(self):
		pass

	@abstractmethod
	def optimize(self):
		"""
		Method to optimize the compute process
		Returns the min OF value in float format and the respective cells in np.array format
		"""
		pass

class Genetic_Algorithm(Optimizer):

	def optimize(self, function, min_value, max_value, pop_size, max_iter, args= [0.05, 0.1]):
		"""
		Method to optimize the compute process using the genetic algorithm

		Returns the min OF value in float format and the respective cells in np.array format
		"""
		mutation = args[0]
		elit = args[1]
		varbound = np.vstack((min_value, max_value)).T
		# varbound = np.array([[0,well.nb_data]]*(simwell.nb_cells-1))
		algorithm_param = {'max_num_iteration': None,\
					   'population_size':pop_size,\
					   'mutation_probability': mutation,\
					   'elit_ratio': elit,\
					   'crossover_probability': 0.5,\
					   'parents_portion': 0.3,\
					   'crossover_type':'uniform',\
					   'max_iteration_without_improv': max_iter}

		model = ga(function=function,dimension=varbound.shape[0],variable_type='real',variable_boundaries=varbound,algorithm_parameters=algorithm_param)
		model.run()

		return model.best_variable, model.best_function

# class Particle_Swarm_Optimization(Optimizer):
#
# 	def optimize(self, function, min_value, max_value, pop_size, max_iter, args=None):
# 		"""
# 		Method to optimize the compute process using the particle swarm optimization
#
# 		Returns the min OF value in float format and the respective cells in np.array format
# 		"""
# 		best_variable, OF = pso(function, min_value, max_value, ieqcons=[], f_ieqcons=None,
# 								args=(), kwargs={},
# 								swarmsize=pop_size, omega=0.7, phip=0.55, phig=0.65, maxiter=max_iter, minstep=1e-7,
# 								minfunc=1e-8, debug=False)
# 		return best_variable, OF


class Particle_Swarm_Optimization(Optimizer):

	def optimize(self, function, min_value, max_value, pop_size, max_iter, args=None):
		"""
		Method to optimize the compute process using the particle swarm optimization

		Returns the min OF value in float format and the respective cells in np.array format
		"""
		omegavalue = args[0]
		phipvalue = args[1]
		phigvalue = args[2]
		best_variable, OF = pso(function, min_value, max_value, ieqcons=[], f_ieqcons=None,
								args=(), kwargs={},
								swarmsize=pop_size, omega=omegavalue, phip=phipvalue, phig=phigvalue, maxiter=max_iter, minstep=1e-7,
								minfunc=1e-8, debug=False)
		return best_variable, OF


class SCOOF_Brute_Force(Optimizer):

	def optimize(self, function, max_value, min_value, pop_size, max_iter, args=None):
		"""
		Method to compute the value of the OF for a geological well

		Returns the min OF value in float format and the respective cells
		"""
		if args == None:
			sys.exit('Error: SCOOFBF can only be used as a internal optimizer!')
		well = args.well
		simwell = args.simwell
		alpha = args.alpha
		beta = args.beta

		# Initialization
		min_OF = 1E8
		interp_well = []

		facies_simwell = np.array(simwell.normdata[well.facies_list])
		simthick = np.array(simwell.normdata['thickness'])

		OF_Table = {}

		for simcell in range(simwell.nb_cells):
			OF = [np.multiply(well.geocells_matrix[i], facies_simwell[simcell]).sum(axis=1) * alpha
				  + abs(well.geothick[i] - simthick[simcell]) * beta for i in range(well.nb_data + 1)]
			for i, ligne in enumerate(OF):
				for j, col in enumerate(ligne):
					OF_Table[(simcell, i, j)] = col

		positions = range(well.nb_data + 1)
		list_cells = range(simwell.nb_cells)
		last = well.nb_data

		# Evaluation of all combinations
		for cells in itertools.combinations_with_replacement(positions, simwell.nb_cells - 1):
			cells = (0,) + cells + (last,)
			OF = 0
			for i in list_cells:
				OF += OF_Table[(i, cells[i], cells[i + 1] - cells[i])]
				if OF >= min_OF:
					break

			if min_OF > OF:
				min_OF = OF
				my_cells = cells

		well_thick = [well.geothick[my_cells[i]][my_cells[i + 1] - my_cells[i]] for i in range(simwell.nb_cells)]
		interpreted_well = [well.geocells_matrix[my_cells[i]][my_cells[i + 1] - my_cells[i]] for i in
							range(simwell.nb_cells)]
		interp_well = interpreted_well * np.array(well_thick)[:, None]

		return np.asarray(interp_well), min_OF

class Dijkstra(Optimizer):

	def optimize(self, function, max_value, min_value, pop_size, max_iter, args=None):
		"""
		Method to optimize the compute process using dijkstra

		Returns the min OF value in float format and the respective cells in np.array format
		"""
		if args == None:
			sys.exit('Error: Dijkstra can only be used as a internal optimizer!')
		well = args.well
		simwell = args.simwell
		alpha = args.alpha
		beta = args.beta

		facies_simwell = args.facies_simwell
		simthick = args.simthick
		wnbd = well.nb_data
		snbc = simwell.nb_cells

		# Construction of a sparse matrix containing the distances between the nodes of the graph
		# Fill with first cell (always starts at first position)
		values = np.multiply(well.geocells_matrix[0], facies_simwell[0]).sum(axis=1) * alpha + abs(
			well.geothick[0] - simthick[0]) * beta
		values = list(values)
		row_ind = [0] * (wnbd + 1)
		col_ind = [*range(1, (wnbd + 1) + 1)]

		# TO DO: mettre geocells_matrix and geothick in a crs_matrix ou brs_matrix pour
		# 	1- gagner du temps au calcul de la FO (multiplication de matrices)
		#	2- ne pas avoir a transformer des listes en brs_matrix ou brs_matrix

		# Fill intermediate cells
		for simcell in range(snbc - 2):
			# Compute the OF values for each dionisos with each possible well data sequence
			for i in range(wnbd + 1):
				values += list(
					np.multiply(well.geocells_matrix[i], facies_simwell[simcell + 1]).sum(axis=1) * alpha + abs(
						well.geothick[i] - simthick[simcell + 1]) * beta)
				row_ind += [simcell * (wnbd + 1) + i + 1] * ((wnbd + 1) - i)
				add_col = [*range((simcell + 1) * (wnbd + 1) + 1 + i, (simcell + 2) * (wnbd + 1) + 1)]
				col_ind += add_col

		# Fill with last cell (always ends at last position)
		values += list(np.multiply(well.geocells_matrix[i][-1], facies_simwell[-1]).sum() * alpha + abs(
			well.geothick[i][-1] - simthick[-1]) * beta for i in range(wnbd + 1))
		row_ind += [2 + (snbc - 1) * (wnbd + 1) - (wnbd - i + 1 + 1) for i in range(wnbd + 1)]
		col_ind += [2 + (snbc - 1) * (wnbd + 1) - 1] * ((wnbd + 1))

		# adding a last point to get a square matrix
		values += [0]
		row_ind += [2 + (snbc - 1) * (wnbd + 1) - 1]
		col_ind += [2 + (snbc - 1) * (wnbd + 1) - 1]
		graph = coo_matrix((values, (row_ind, col_ind)))

		# Computing the shortest path
		dist_matrix, predecessors = dijkstra(csgraph=graph, directed=True, indices=0, return_predecessors=True)

		# reconstruction of the shortest path (intepreted well)
		markers = [(snbc - 1) * (wnbd + 1) + 1]
		for i in range(snbc):
			markers = [predecessors[markers[-1 - i]]] + markers
		markers[1:-1] = [(i - 1) % (wnbd + 1) for i in markers[1:-1]]
		markers[-1] = wnbd

		well_thick = [well.geothick[markers[i]][markers[i + 1] - markers[i]] for i in range(simwell.nb_cells)]
		interpreted_well = [well.geocells_matrix[markers[i]][markers[i + 1] - markers[i]] for i in
							range(simwell.nb_cells)]
		interp_well = interpreted_well * np.array(well_thick)[:, None]

		return np.asarray(markers, dtype=np.float32), dist_matrix[-1]

class Dijkstra_With_Square_Root_Thickness(Optimizer):

	def optimize(self, function, max_value, min_value, pop_size, max_iter, args=None):
		"""
		Method to optimize the compute process using dijkstra

		Returns the min OF value in float format and the respective cells in np.array format
		"""

		if args == None:
			sys.exit('Error: Dijkstra can only be used as a internal optimizer!')
		well = args.well
		simwell = args.simwell
		alpha = args.alpha
		beta = args.beta



		facies_simwell = args.facies_simwell
		simthick = args.simthick
		wnbd = well.nb_data
		snbc = simwell.nb_cells

		# well_geothick_squareroot = [thickness_value ** 2 for thickness_value in well.geothick[0]]
		# simthick_squareroot = [thickness_value ** 2 for thickness_value in simthick[0]]

		# Construction of a sparse matrix containing the distances between the nodes of the graph
		# Fill with first cell (always starts at first position)
		values = np.multiply(well.geocells_matrix[0], facies_simwell[0]).sum(axis=1) * alpha + abs(([i ** (1/2) for i in well.geothick[0]])
			 - (simthick[0]**(1/2)))* beta
		values = list(values)
		row_ind = [0] * (wnbd + 1)
		col_ind = [*range(1, (wnbd + 1) + 1)]

		# TO DO: mettre geocells_matrix and geothick in a crs_matrix ou brs_matrix pour
		# 	1- gagner du temps au calcul de la FO (multiplication de matrices)
		#	2- ne pas avoir a transformer des listes en brs_matrix ou brs_matrix

		# Fill intermediate cells
		for simcell in range(snbc - 2):
			# Compute the OF values for each dionisos with each possible well data sequence
			for i in range(wnbd + 1):
				values += list(
					np.multiply(well.geocells_matrix[i], facies_simwell[simcell + 1]).sum(axis=1) * alpha + abs(([i ** (1/2) for i in well.geothick[i]])
						 - (simthick[simcell + 1]**(1/2))) * beta)
				row_ind += [simcell * (wnbd + 1) + i + 1] * ((wnbd + 1) - i)
				add_col = [*range((simcell + 1) * (wnbd + 1) + 1 + i, (simcell + 2) * (wnbd + 1) + 1)]
				col_ind += add_col

		# Fill with last cell (always ends at last position)
		values += list(np.multiply(well.geocells_matrix[i][-1], facies_simwell[-1]).sum() * alpha + abs((well.geothick[i][-1]**(1/2))
			 - (simthick[-1]**(1/2))) * beta for i in range(wnbd + 1))
		row_ind += [2 + (snbc - 1) * (wnbd + 1) - (wnbd - i + 1 + 1) for i in range(wnbd + 1)]
		col_ind += [2 + (snbc - 1) * (wnbd + 1) - 1] * ((wnbd + 1))

		# adding a last point to get a square matrix
		values += [0]
		row_ind += [2 + (snbc - 1) * (wnbd + 1) - 1]
		col_ind += [2 + (snbc - 1) * (wnbd + 1) - 1]
		graph = coo_matrix((values, (row_ind, col_ind)))

		# Computing the shortest path
		dist_matrix, predecessors = dijkstra(csgraph=graph, directed=True, indices=0, return_predecessors=True)

		# reconstruction of the shortest path (intepreted well)
		markers = [(snbc - 1) * (wnbd + 1) + 1]
		for i in range(snbc):
			markers = [predecessors[markers[-1 - i]]] + markers
		markers[1:-1] = [(i - 1) % (wnbd + 1) for i in markers[1:-1]]
		markers[-1] = wnbd

		well_thick = [well.geothick[markers[i]][markers[i + 1] - markers[i]] for i in range(simwell.nb_cells)]
		interpreted_well = [well.geocells_matrix[markers[i]][markers[i + 1] - markers[i]] for i in
							range(simwell.nb_cells)]
		interp_well = interpreted_well * np.array(well_thick)[:, None]

		return np.asarray(markers, dtype=np.float32), dist_matrix[-1]

class Brute_Force(Optimizer):

	def optimize(self, function, min_value, max_value, pop_size, max_iter, args=None):
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

		if min_value.size != max_value.size:
			sys.exit("Error: Lower and upper bounds sizes must be the same!")

		if max_iter < 2:
			sys.exit("Error: Number of samples must be equal or greater than 2!")

		uncertain_parameters_sample = []
		for parameter in range(min_value.size):
			uncertain_parameters_sample.append(np.linspace(min_value[parameter], max_value[parameter], max_iter))

		uncertain_parameters_combinations = itertools.product(*tuple(uncertain_parameters_sample))

		total_OF_Values = []
		for combination in uncertain_parameters_combinations:
			total_OF_Values.append(function(np.asarray(combination)))

		uncertain_parameters_list = list(itertools.product(*tuple(uncertain_parameters_sample)))

		best_variable = uncertain_parameters_list[total_OF_Values.index(min(total_OF_Values))],
		OF = min(total_OF_Values)

		return best_variable, OF