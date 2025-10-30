"""
Main file of the prototype of calculation of calibration workflow -
Calibration of DionisosFlow models
"""
# Import libraries
import pandas as pd
from sklearn.model_selection import ParameterGrid, ParameterSampler
from scipy.stats import uniform
import os
import time
import shutil
import Code_CW.Objects as Objects
import Code_CW.Functions as Functions

# NOVAS IMPORTAÇÕES
from skopt import gp_minimize
from skopt.space import Real
from skopt.utils import use_named_args
import matplotlib.pyplot as plt # Para o gráfico de convergência

#########################
tic = time.time()
#########################
#       INPUT DATA     #
#########################
# OpenFlow Dionisos Simulator (.exe) file needed to run a dionisos simulation.
#enginefolder = '"C:\\Program Files\\Beicip\\OpenFlowSuite_2021\\WfEngine\plugins\\fr.ifp.dionisosflow.arcades.engine.win64_13.1.0.r25519\\windows\\ArcaDESLauncher.exe"'
enginefolder = r'"C:\Program Files\Beicip\OpenFlowSuite_2023\WfEngine\plugins\fr.ifp.dionisosflow.arcades.engine.win64_15.1.0.r27540\windows\ArcaDESLauncher.exe"'
## FILES AND FOLDERS
# Main project directory/working Folder
# ROOT BOING
#
#root_folder = "C:\\Users\\09008889978\\Desktop\\WorkingDir"
# ROOT LINDAURA
#root_folder = "C:\\Users\\00615045995\\PycharmProjects\\Calibration_workflow"
# ROOT JESSICA
root_folder = r"C:\Users\00584647271\Desktop\diretorio"
#root_folder = "C:\\Users\\45597483811\\PycharmProjects\\Calibration_workflow_v2.0_refatorado"

### STUDY CASE
case = "Hypertuning_Bayesian_optimized"

### DIONISOS PROJECT INPUT
# Absolute initial coordinates of the project's grid defined on Dionisos interface.
# x0 = 405000
# y0 = 1055000
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
#Mixed Case Lithologies
litho_class_ids = {'1':'Silt','2':'Sand','3':'Reef','4':'Carbonate'}
carbo_ids = ['3','4']
mud_text_ids = ['1']
sand_text_ids = ['2','3']
gravel_text_ids = ['4']

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
facies_def_file = os.path.join(simdir, 'Table_def_facies.txt')
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
internal_optimization_method = 'DSR' #'PSO'= Particle Swarm Optimization or 'GA'= Genetic Algorithm or 'SCOOFBF' = SCOOF with brute force,'D'= Dijkstra or 'DSR'= Dijkstra with with square root in the thicknesses
external_optimization_method = 'PSO' #'PSO'= Particle Swarm Optimization or 'AG'= Genetic Algorithm or 'BF'= Brute Force
alpha = 4
beta = 1
                                        #PROOF PARAMETER TO SET REFERENCE SOLUTION

pop_size = 1 # Population size
max_iter = 1 # Maximum number of iterations (termination criterion)


# ==========================================================
""" DEFINA SEU MÉTODO DE CALIBRAÇÃO AQUI """
# REMOVA O BLOCO DE GRID/RANDOM SEARCH
# if CALIBRATION_METHOD == 'G':
#     ...
# elif CALIBRATION_METHOD == 'R':
#     ...
# ==========================================================

# 1. DEFINA O ESPAÇO DE BUSCA PARA A OTIMIZAÇÃO BAYESIANA
search_space = [
    Real(0.0, 1.0, name='omega'),
    Real(0.0, 1.0, name='phip'),
    Real(0.0, 1.0, name='phig'),
]

# 2. CRIE UMA LISTA PARA ARMAZENAR OS RESULTADOS DE CADA CHAMADA
results_list = []
# Contador global para nomear as pastas de resultado de forma única
run_counter = 0


# 3. CRIE A FUNÇÃO OBJETIVO QUE SERÁ CHAMADA PELO gp_minimize

@use_named_args(search_space)
def objective_function(omega, phip, phig):
    """
    Esta função executa UMA simulação completa com um conjunto de hiperparâmetros
    e retorna o valor da função objetivo (OF_value) para ser minimizado.
    """
    global run_counter
    run_counter += 1

    print(f"\n--- Iniciando Iteração Bayesiana Nº {run_counter} ---")
    print(f"Parâmetros: omega={omega:.3f}, phip={phip:.3f}, phig={phig:.3f}")

    # Definição da pasta de destino para esta combinação
    results_folder_final = os.path.join(simdir, f"Results_{run_counter}")
    default_results_folder = os.path.join(simdir, "Results")

    # Limpeza de pastas
    if os.path.exists(results_folder_final):
        shutil.rmtree(results_folder_final)
    if os.path.exists(default_results_folder):
        shutil.rmtree(default_results_folder)

    # A lógica abaixo é a mesma do seu loop 'for' original
    simulation = Objects.Simulation_parameters()
    internal_optimizer = Functions.get_optimizer(internal_optimization_method)
    external_optimizer = Functions.get_optimizer(external_optimization_method)

    OF = Functions.get_objective_function(OF_type)
    OF.set_objective_function(case, simdir, arcfilebase, origin, well_files, extract_columns, facies_definition,
                              simulation,
                              litho_class_ids=litho_class_ids, carbo_ids=carbo_ids, mud_ids=mud_text_ids,
                              sand_ids=sand_text_ids, gravel_ids=gravel_text_ids,
                              alpha=alpha, beta=beta,
                              internal_optimization_method=internal_optimization_method)
    OF.set_optmization(uncertain_parameters, enginefolder, internal_optimizer, external_optimizer,
                       reference_value=[])

    tic_run = time.time()
    best_variable, OF_value = external_optimizer.optimize(
        OF.compute,
        uncertain_parameters['min_value'].to_numpy(),
        uncertain_parameters['max_value'].to_numpy(),
        pop_size,
        max_iter,
        [omega, phip, phig]
    )
    tac_run = time.time() - tic_run

    # Processamento de resultados (igual ao seu código)
    Functions.plotResults(OF_type, simulation, simdir)
    # Limpeza de arquivos Excel
    # if os.path.exists(default_results_folder):
    #     for root, dirs, files in os.walk(default_results_folder, topdown=False):
    #         for name in files:
    #             if name.endswith(('.xlsx', '.xls')):
    #                 os.remove(os.path.join(root, name))
    # Mover pasta de resultados
    if os.path.exists(default_results_folder):
        shutil.move(default_results_folder, results_folder_final)

    # Armazenar resultados desta iteração
    current_result = {
        'Combination': run_counter,
        'omega': omega,
        'phip': phip,
        'phig': phig,
        'OF_value': OF_value,
        'time_h': tac_run / 3600
    }
    for param_name, param_value in zip(uncertain_parameters["uncertain_parameter"], best_variable):
        current_result[param_name] = param_value

    results_list.append(current_result)

    # # Salva um backup a cada iteração
    # df_results = pd.DataFrame(results_list)
    # output_file = os.path.join(root_folder, case, "calibration_results_bayesian.xlsx")
    # df_results.to_excel(output_file, index=False)

    print(f"Tempo da iteração: {tac_run:.2f} s | OF_value: {OF_value:.4f}")
    print(f"Melhor conjunto de parâmetros encontrado: {best_variable}")
    print(f"--- Fim da Iteração Nº {run_counter} ---")

    # A função DEVE retornar o valor a ser minimizado
    return OF_value

# =====================================================================================================================
"""REMOVA TODO O SEU ANTIGO LOOP DE SIMULAÇÃO "for idx, row in df.iterrows():"""
# =====================================================================================================================


# 4. EXECUTE A OTIMIZAÇÃO BAYESIANA
N_CALLS = 10 # Número total de simulações a serem executadas (ajuste conforme necessário)
print(f"\nIniciando Otimização Bayesiana com {N_CALLS} chamadas...")

result_bayesian = gp_minimize(
    func=objective_function,
    dimensions=search_space,
    n_calls=N_CALLS,
    verbose=True
)

# 5. PROCESSAR E EXIBIR OS RESULTADOS FINAIS
print("\n--- Otimização Bayesiana Concluída ---")
print(f"Melhor valor de OF encontrado: {result_bayesian.fun:.4f}")
print("Melhores hiperparâmetros (omega, phip, phig):")
best_hyperparams = {name: val for name, val in zip(['omega', 'phip', 'phig'], result_bayesian.x)}
print(best_hyperparams)

# Salvar o DataFrame final com todos os resultados
output_file = os.path.join(root_folder, case, "calibration_results_bayesian.xlsx")
final_df = pd.DataFrame(results_list)
final_df.to_excel(output_file, index=False)
print(f"\nResultados completos salvos em: {output_file}")

# Opcional: Plotar a convergência
from skopt.plots import plot_convergence

plot_convergence(result_bayesian)
convergence_plot_path = os.path.join(root_folder, case, "convergence_plot.png")
plt.savefig(convergence_plot_path)
# print(f"Gráfico de convergência salvo em: {convergence_plot_path}")
# plt.show()




























