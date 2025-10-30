import time
import pandas as pd
from sklearn.model_selection import ParameterGrid, ParameterSampler
from scipy.stats import uniform
import BayesianOptimization

""""" Escolha o método: GRID ou RANDOM """""

# ## --- GRID SEARCH ---
# param_grid = {
#     'omega': [0, 0.25, 0.5, 0.75, 1],
#     'phip': [0, 0.25, 0.5, 0.75, 1],
#     'phig': [0, 0.25, 0.5, 0.75, 1],
# }
# df = pd.DataFrame(list(ParameterGrid(param_grid)))
#combinations = df.to_dict(orient="records")   # transforma em lista de dicts

## --- RANDOM SEARCH ---
param_dist = {
    'omega': uniform(0, 1),
    'phip': uniform(0, 1),
    'phig': uniform(0, 1),
}
random_combinations = list(ParameterSampler(param_dist, n_iter=10, random_state=42))
combinations = [{k: round(v, 2) for k, v in params.items()} for params in random_combinations]
df = pd.DataFrame(combinations)

### LOOP DE SIMULAÇÃO ###
# for params in combinations:
#     omega = params['omega']
#     phip = params['phip']
#     phig = params['phig']

#######

for idx, row in df.iterrows():
    omega = row['omega']
    phip = row['phip']
    phig = row['phig']

# row = df.iloc[4]
# print(row)