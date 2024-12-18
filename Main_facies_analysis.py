import os
import lasio
import time
import Code_CW.Analysis as Analysis

tic = time.time()
#Root folder containing necessary files and folders
root = 'C:\\Users\\joaov\\Desktop\\Artigo\\Projeto\\facies_analysis'

#Folder containing well files
wells_dir = os.path.join(root, 'wells')

#File containing markers of the top and bottom of the zone of interest
markers_file = os.path.join(root,'wellMarkers.txt')

color_ref_file = os.path.join(root, 'color_reference.txt')
#Data columns to be read
depth_column = 'DEPT'
facies_column = 'fac'
bathymetry_column = 'bat'
litho_column = 'lito_upscaled'
extract_columns = [depth_column, facies_column, bathymetry_column, litho_column]
del facies_column, depth_column, bathymetry_column, litho_column

#Domain [x0, xlen, y0, ylen]
domain = [808148, 8400, 8916612, 7500]

#Analysis Functions
FaciesAnalysis = Analysis.FaciesAnalysis(root, wells_dir, markers_file, extract_columns, domain, color_ref_file)
FaciesAnalysis.wellPieCharts(FaciesAnalysis.totalinfo, FaciesAnalysis.faciescolor)
FaciesAnalysis.scatterpie(FaciesAnalysis.totalinfo)

tac = time.time()
time = tac - tic
print(f'Facies Analysis Ended. Computational time: {time}s.')