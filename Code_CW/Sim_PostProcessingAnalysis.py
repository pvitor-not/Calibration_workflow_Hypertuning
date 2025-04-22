import pandas as pd
import Code_CW.Functions as Functions
from pandas import DataFrame as df
import os
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.patches import Patch
import warnings
from fpdf import FPDF
from datetime import datetime
import lasio
import numpy as np
import Code_CW.StatisticalAnalysis as Statistics


class PropSimAnalysis:
    def __init__(self, root, of_type, well_files, id_file, extract_columns, domain,
                 color_reference_file, facies_color):
        warnings.simplefilter(action='ignore', category=FutureWarning)
        self.root = root
        self.of_type = of_type.upper()
        self.id_file = id_file
        self.id_table = self.readIDFile(id_file)
        self.extract_columns = extract_columns
        self.domain = domain
        self.setdomain()
        self.sedProp = self.set_color_map(color_reference_file)
        self.facies_color = self.set_colorfacies_map(facies_color)
        self.prop_columns = self.set_proportion_columns()
        self.cell_reference = self.getCellReference()
        self.wells_file = well_files
        self.root_result = self.getRootResultArq()
        self.welldataprop = self.getWellData(self.wells_file, self.extract_columns, self.prop_columns, self.root_result,
                                             self.cell_reference, self.of_type)

    def getCellReference(self):
        if self.of_type == "S":
            return "SimWell"
        elif self.of_type == "P":
            return ""
        else:
            BaseException("Error! OF_type must be set as 'S' or 'P.'")
            raise "Error! OF_type must be set as 'S' or 'P.'"

    def setdomain(self):
        self.xo = self.domain[0]
        self.xlen = self.domain[1]
        self.yo = self.domain[2]
        self.ylen = self.domain[3]

    def readIDFile(self, id_file):
        if id_file != '':
            f = open(id_file, 'r')
            # getting the number of properties
            data = ['SubZone', 'Zone']

            if f.readline()[:-1] != "ID":
                print('Error: the file is not a file of definition of zone ID ')
            else:
                try:
                    nb_id = int(f.readline())
                except Exception as e:
                    print('Exception:', e)
            # Getting the list of properties
            for id in range(nb_id):
                data.append(f.readline()[:-1])

            # Getting the number of zones
            if f.readline()[:-1] != "Zone":
                print('Error: You must set a section for definition of number of zones below the ID section')
            else:
                try:
                    nb_zones = int(f.readline())
                except Exception as e:
                    print('Exception:', e)

                zone_name = f.readline()[:-1]

            # getting the number of subzones
            if f.readline()[:-1] != "SubZone":
                print('Error: You must set a section for definition of number of subzones below the zone section')
            else:
                try:
                    nb_subzones = int(f.readline())
                except Exception as e:
                    print('Exception:', e)
            if nb_subzones == 0:
                data_to_table = [[0, zone_name, '-', '-']]
                id_table = df(data=data_to_table, columns=data)

                f.close()
            else:

                f.close()

                id_table = pd.read_csv(id_file, delim_whitespace=True, header=None, skiprows=(nb_id + 7))
                id_table.columns = data
        else:
            id_table = ''
        return id_table

    def set_color_map(self, color_reference_file):
        if color_reference_file != '':
            colors_ref = pd.read_csv(color_reference_file, delim_whitespace=True, header=0)
            # colors_ref.sort_values(by='Facies', inplace=True)
            for i in range(0, len(colors_ref.SedProp)):
                colors_ref.SedProp.loc[i] = str(colors_ref.SedProp.loc[i])
            colors_ref.set_index('SedProp', inplace=True)
            colors_ref.dropna(inplace=True)
        else:
            colors_ref = ''
        return colors_ref

    def set_colorfacies_map(self, color_reference_file):
        if color_reference_file != '':
            colors_ref = pd.read_csv(color_reference_file, delim_whitespace=True, header=0, skiprows=1)
            for i in range(0, len(colors_ref.Facie)):
                colors_ref.Facie.loc[i] = str(colors_ref.Facie.loc[i])
            colors_ref.set_index('Facie', inplace=True)
            colors_ref.dropna(inplace=True)
        else:
            colors_ref = ''
        return colors_ref

    def set_proportion_columns(self):
        column_names = [x for x in self.sedProp.index]
        return column_names

    def getRootResultArq(self):
        of_results = os.path.join(f'{self.root}\\Results\\{self.of_type}\\OFresults.xlsx')
        try:
            nb_BestSim = pd.read_excel(of_results, 'OF_values', header=None)[0][1].split(" ")
        except:
            raise "Error! Variable 'OF_type' is set incorrect."

        if self.of_type == "S":
            simWellPath = os.path.join(
                f'{self.root}\\Results\\{self.of_type}\\Sim_{nb_BestSim[1]}_Well_Facies_Sequences.xlsx')
        else:
            simWellPath = os.path.join(f'{self.root}\\Results\\{self.of_type}\\PROOF_Stats\\Sim{nb_BestSim[1]}'
                                       f'\\Well_Sequences.xlsx')
        return simWellPath

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

    def simPieCharts(self):
        if self.of_type == "S":
            return False

    def getWellData(self, files, extract_columns, prop_columns, root_result, cell_reference, of_type):
        welldata = {}
        for well in files:
            wellinfo = Functions.readWellFile(well, extract_columns)
            welldata[f"{wellinfo['name']}"] = PAGeoWell(wellinfo, prop_columns, root_result, cell_reference, of_type)
        return welldata

    def normalizeData(self, dataDict):
        soma = 0
        for value in dataDict.items():
            soma += value[1]

        for key, item in dataDict.items():
            dataDict[key] = round(dataDict[key] / soma, ndigits=4)

        keys_to_extract = [key for key, item in dataDict.items() if dataDict[key] <= 0.009]

        for i in keys_to_extract:
            del dataDict[i]

        return dataDict

    def plot_wellpiechart(self, values, labels, colors, well, subzone):
        def formatpctvalues(data):
            pct = '{:.1f}%'.format(data)
            return pct

        if subzone is not None:  # If you do calibrate using subzones!!!
            # Single well-PieChart for each subzone with the label and values.
            plt.pie(values, labels=labels, colors=colors, pctdistance=0.65,
                    autopct=formatpctvalues,
                    startangle=0, radius=.75,
                    textprops=dict(size=8, fontweight='bold'))
            plt.title(f"Well {well} - Subzone {subzone}")
            plt.savefig(f"{self.root}\\Results\\{self.of_type}\\WellPieChart\\Single_WellPieCharts\\Well {well}"
                        f" - Subzone {subzone}.png", transparent=True)
            plt.close()

            # Single well-PieChart for each subzone without the label and values and bigger for scatter in map.
            plt.pie(values, labels=None, colors=colors, pctdistance=0.65,
                    startangle=0, radius=1.5)
            plt.title(f"{well}", pad=1, y=.5, fontdict=dict(size=55, fontweight='bold'))
            plt.savefig(f"{self.root}\\Results\\{self.of_type}\\WellPieChart\\PieCharts_toMap\\Well {well} "
                        f"Subzone {subzone}.png", transparent=True)
            plt.close()
        else:  # If you do not calibrate using subzones!!!
            # Single well-PieChart for the entire well section with the label and values.
            plt.pie(values, labels=labels, colors=colors, pctdistance=0.65,
                    autopct=formatpctvalues,
                    startangle=0, radius=.75,
                    textprops=dict(size=8, fontweight='bold'))
            plt.title(f"Well {well} - Entire Well Section")
            plt.savefig(
                f"{self.root}\\Results\\{self.of_type}\\WellPieChart\\Single_WellPieCharts\\Well {well} - "
                f"Entire Well Section.png", transparent=True)
            plt.close()

            # Single well-PieChart for the entire section without the label and values and bigger for scatter in map.
            plt.pie(values, labels=None, colors=colors, pctdistance=0.65,
                    startangle=0, radius=1.5)
            plt.title(f"{well}", pad=1, y=.5, fontdict=dict(size=55, fontweight='bold'))
            plt.savefig(
                f"{self.root}\\Results\\{self.of_type}\\WellPieChart\\PieCharts_toMap\\Well {well} - "
                f"Entire Well Section.png", transparent=True)
            plt.close()

    def wellPieCharts(self, wellDF, id_table, color_reference):
        if not os.path.isdir(f'{self.root}\\Results\\{self.of_type}\\WellPieChart'):
            os.makedirs(f'{self.root}\\Results\\{self.of_type}\\WellPieChart')
        if not os.path.isdir(f'{self.root}\\Results\\{self.of_type}\\WellPieChart\\Single_WellPieCharts'):
            os.makedirs(f'{self.root}\\Results\\{self.of_type}\\WellPieChart\\Single_WellPieCharts')
        if not os.path.isdir(f'{self.root}\\Results\\{self.of_type}\\WellPieChart\\PieCharts_toMap'):
            os.makedirs(f'{self.root}\\Results\\{self.of_type}\\WellPieChart\\PieCharts_toMap')

        for well, data in wellDF.items():

            if len(id_table.index) == 1:  # If you do not calibrate using subzones!!!
                data_to_plot = {}
                for key in data.dfProportions.keys():
                    data_to_plot[key] = data.dfProportions[key].sum()

                data_to_plot = self.normalizeData(data_to_plot)
                wellcolors = [color_reference.loc[x] for x in color_reference.index if x in data_to_plot.keys()]
                self.plot_wellpiechart(data_to_plot.values(), data_to_plot.keys(), wellcolors, well, subzone=None)

            else:  # If you do calibrate using subzones!!!
                # plot for the subzones well section
                for idx in id_table.index:
                    data_to_plot = {}

                    for key in data.dfProportions.keys():
                        data_to_plot[key] = data.dfProportions[key][
                                            id_table['Start'][idx]:id_table['Finish'][idx] + 1].sum()

                    data_to_plot = self.normalizeData(data_to_plot)

                    wellcolors = [color_reference.loc[x] for x in color_reference.index if x in data_to_plot.keys()]
                    self.plot_wellpiechart(data_to_plot.values(), data_to_plot.keys(), wellcolors, well,
                                           subzone=id_table['SubZone'][idx])

                # Plot for the entire well section
                data_to_plot = {}
                for key in data.dfProportions.keys():
                    data_to_plot[key] = data.dfProportions[key].sum()

                data_to_plot = self.normalizeData(data_to_plot)
                wellcolors = [color_reference.loc[x] for x in color_reference.index if x in data_to_plot.keys()]
                self.plot_wellpiechart(data_to_plot.values(), data_to_plot.keys(), wellcolors, well, subzone=None)

    def getImage(self, path):
        return OffsetImage(plt.imread(path), zoom=0.12)

    def scatter_for_map(self, data, labels, pietype, subzone):
        fig, ax = plt.subplots(dpi=300, figsize=(8, 6))
        ax.scatter(data.X, data.Y, s=0.02)
        for i in data.index:
            ab = AnnotationBbox(self.getImage(data.PATH.loc[i]), (data.X.loc[i], data.Y.loc[i]),
                                frameon=False, pad=1.3)
            ax.add_artist(ab)
        plt.xlabel('East')
        plt.xlim(self.xo, self.xo + self.xlen)
        plt.ylabel('North')
        plt.ylim(self.yo, self.yo + self.ylen)
        plt.ticklabel_format(axis='y', style='plain')
        plt.legend(handles=labels, fontsize=8)
        if pietype == 'Subzone':
            plt.title(f"Well-PieCharts for Subzone {subzone}")
            plt.savefig(f"{self.root}\\Results\\{self.of_type}\\WellPieChart\\Maps\\_ScatterPie_Map SubZone_"
                        f"{subzone}.png")
        else:
            plt.title(f"Well-PieCharts for the entire well section")
            plt.savefig(f"{self.root}\\Results\\{self.of_type}\\WellPieChart\\Maps\\_ScatterPie_Map _"
                        f"EntireWellSection.png")
        plt.clf()

    def scatterpie(self, wellData, id_table, color_reference):
        path = f'{self.root}\\Results\\{self.of_type}\\WellPieChart\\PieCharts_toMap'
        if not os.path.isdir(f'{self.root}\\Results\\{self.of_type}\\WellPieChart\\Maps'):
            os.makedirs(f'{self.root}\\Results\\{self.of_type}\\WellPieChart\\Maps')

        wellcolors = [(color_reference.loc[x, 'R'], color_reference.loc[x, 'G'], color_reference.loc[x, 'B'])
                      for x in color_reference.index]
        names = [x for x in color_reference.index]
        patches = [Patch(color=c, label=l) for c, l in zip(wellcolors, names)]

        if len(id_table.index) == 1:  # If you do not calibrate using subzones!!!
            wellcoord = df(data=None, columns=['X', 'Y', 'PATH'])
            pies = []
            for r, d, f in os.walk(path):
                for file in f:
                    pies.append(os.path.join(r, file))

            for pie in pies:
                base = os.path.basename(pie)
                wellpie = os.path.splitext(base)
                wellpie = f'{wellpie[0]}'

                well_name = f"{wellpie.split(' ')[1]}"

                wellcoord.loc[f'{format(wellpie)}', 'X'] = float(wellData[well_name].x)
                wellcoord.loc[f'{format(wellpie)}', 'Y'] = float(wellData[well_name].y)
                wellcoord.loc[f'{format(wellpie)}', 'PATH'] = pie

            self.scatter_for_map(data=wellcoord, labels=patches, pietype='AllSection', subzone=None)

        else:  # If you do calibrate using subzones!!!
            # Scatter pie for the subzones!
            for indx in id_table.index:
                wellcoord = df(data=None, columns=['X', 'Y', 'PATH'])
                pies = []
                for r, d, f in os.walk(path):
                    for file in f:
                        file_basename = os.path.splitext(file)
                        file_well_subzone = file_basename[0].split(' ')[-1]
                        if str(file_well_subzone) == f"{id_table['SubZone'][indx]}" \
                                or str(file_well_subzone) == f"{id_table['Zone'][indx]}":
                            pies.append(os.path.join(r, file))

                for pie in pies:
                    base = os.path.basename(pie)
                    wellpie = os.path.splitext(base)
                    wellpie = f'{wellpie[0]}'

                    well_name = f"{wellpie.split(' ')[1]}"

                    wellcoord.loc[f'{format(wellpie)}', 'X'] = float(wellData[well_name].x)
                    wellcoord.loc[f'{format(wellpie)}', 'Y'] = float(wellData[well_name].y)
                    wellcoord.loc[f'{format(wellpie)}', 'PATH'] = pie

                self.scatter_for_map(data=wellcoord, labels=patches, pietype='Subzone',
                                     subzone=id_table['SubZone'][indx])

            # Scatter pie for the entire zone!
            wellcoord = df(data=None, columns=['X', 'Y', 'PATH'])
            pies = []
            for r, d, f in os.walk(path):
                for file in f:
                    file_basename = os.path.splitext(file)
                    file_well_subzone = file_basename[0].split(' ')[-1]
                    if str(file_well_subzone) == f"Section":
                        pies.append(os.path.join(r, file))

            for pie in pies:
                base = os.path.basename(pie)
                wellpie = os.path.splitext(base)
                wellpie = f'{wellpie[0]}'

                well_name = f"{wellpie.split(' ')[1]}"

                wellcoord.loc[f'{format(wellpie)}', 'X'] = float(wellData[well_name].x)
                wellcoord.loc[f'{format(wellpie)}', 'Y'] = float(wellData[well_name].y)
                wellcoord.loc[f'{format(wellpie)}', 'PATH'] = pie

            self.scatter_for_map(data=wellcoord, labels=patches, pietype='AllSection', subzone=None)

    @staticmethod
    def calcular_posicoes(object_class):
        posicoes = [0]
        for i in range(1, len(object_class.df_analysis['Thickness'])):
            posicoes.append(posicoes[-1] + object_class.df_analysis['Thickness'][i - 1])
        object_class.position = posicoes

    def plotar_perfil(self, welldata: dict, cores):
        num_perfis = len(welldata)
        max_depth = max([sum(item.df_analysis['Thickness'][:-1]) for item in welldata.values()]) #########

        fig, axs = plt.subplots(1, num_perfis, figsize=(15, 5), sharey=False)
        fig.suptitle('Perfis de Lithofacies')

        for ax, (key, item) in zip(axs, welldata.items()):
            self.calcular_posicoes(item)

            for i in range(len(item.df_analysis['Thickness'])):
                wellcolors = (cores.loc[f"{item.df_analysis['facies'][i]}", "R"],
                              cores.loc[f"{item.df_analysis['facies'][i]}", "G"],
                              cores.loc[f"{item.df_analysis['facies'][i]}", "B"])

                facies_label = str(item.df_analysis['facies'][i])

                ax.barh(y=item.position[i], width=1, height=item.df_analysis['Thickness'][i], color=wellcolors,
                        align='edge')  # hatch=hatch_pattern)

                ax.text(x=0.5, y=item.position[i] + item.df_analysis['Thickness'][i] / 2,
                        s=facies_label, ha='center', va='center', color='black')

            ax.set_xlabel(key)
            ax.set_ylabel('Profundidade')

            ax.set_ylim([max_depth, 0])

            perfil_max_depth = sum(item.df_analysis['Thickness'])
            step = perfil_max_depth / 10.1

            yticks = np.arange(0, perfil_max_depth, 30)
            if perfil_max_depth not in yticks:
                yticks = np.append(yticks, perfil_max_depth)
            ax.set_yticks(yticks)
            ax.set_xticks([])
            ax.set_xticklabels([])

            for spine in ax.spines.values():
                spine.set_visible(False)

            ax.yaxis.set_tick_params(which='both', labelleft=True)

        path_result = f'{self.root}\\Results\\{self.of_type}'
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        plt.savefig(f'{path_result}\\_Wells_Log.png')
        plt.close()

    def calcular_posicoesOBS(self, welldata):
        sequence_df = welldata.sequence
        tck = []
        fafaci = sequence_df['facies'].tolist()
        thickness = sequence_df['Thickness'].tolist()
        faciesobs = []

        i = 0
        while i < len(fafaci):
            facie_atual = fafaci[i]

            sumtck = 0
            while i < len(fafaci) and fafaci[i] == facie_atual:
                sumtck += thickness[i]
                i += 1
            faciesobs.append(facie_atual)
            tck.append(sumtck)

        data = {
            'facies': faciesobs,
            'Thickness': tck
        }
        posit_df = pd.DataFrame(data)

        posicoes = [0]
        for i in range(1, len(posit_df['Thickness'])):
            posicoes.append(posicoes[-1] + posit_df['Thickness'][i - 1])
        posit_df['position'] = posicoes

        return posit_df

    def plotar_perfilOBS(self, welldata: dict, cores):
        num_perfis = len(welldata)
        max_depth = max([sum(item.sequence['Thickness']) for item in welldata.values()])
        fig, axs = plt.subplots(1, num_perfis, figsize=(15, 5), sharey=False)
        fig.suptitle('Perfis de Lithofacies')

        for ax, (key, item) in zip(axs, welldata.items()):
            posit_df = self.calcular_posicoesOBS(item)

            for i in range(len(posit_df['facies'])):
                wellcolors = (cores.loc[f"{int(posit_df['facies'][i])}", "R"],
                              cores.loc[f"{int(posit_df['facies'][i])}", "G"],
                              cores.loc[f"{int(posit_df['facies'][i])}", "B"])

                facies_label = str(posit_df['facies'][i])

                ax.barh(y=posit_df['position'][i], width=1, height=posit_df['Thickness'][i], color=wellcolors,
                        align='edge')
                ax.text(x=0.5, y=posit_df['position'][i] + posit_df['Thickness'][i] / 2,
                        s=facies_label, ha='center', va='center', color='black')

            ax.set_xlabel(key)
            ax.set_ylabel('Profundidade')

            ax.set_ylim([max_depth, 0])

            perfil_max_depth = sum(posit_df['Thickness'])
            step = perfil_max_depth / 10.1
            yticks = np.arange(0, perfil_max_depth, 30)
            if perfil_max_depth not in yticks:
                yticks = np.append(yticks, perfil_max_depth)
            ax.set_yticks(yticks)
            ax.set_xticks([])
            ax.set_xticklabels([])

            for spine in ax.spines.values():
                spine.set_visible(False)

            ax.yaxis.set_tick_params(which='both', labelleft=True)

        path_result = f'{self.root}\\Results\\{self.of_type}'
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(f'{path_result}\\_WellsOBS_Log.png')
        plt.close()

    def faciesExtraction(self, wells, root, cell_ref, markers):
        wells_perfil = {}
        marker_table = self.readMarkersFile(markers)
        for well in wells:
            file = lasio.read(well)
            wells_perfil[file.well.WELL.value] = PerfilFacies(file, marker_table, cell_ref, root)

        self.plotar_perfil(wells_perfil, self.facies_color)
        self.plotar_perfilOBS(wells_perfil, self.facies_color)


class PAGeoWell:
    def __init__(self, wellinfo, prop_columns, root_result, cell_reference, of_type):
        self.name = wellinfo['name']
        self.x = wellinfo['x']
        self.y = wellinfo['y']
        self.dfProportions = self.getSedProportion(root_result, cell_reference, prop_columns, of_type)

    def getSedProportion(self, root, cell_reference, prop_columns, of_type):
        f = pd.read_excel(root, sheet_name=f'{cell_reference}{self.name}')
        columns = f.columns
        dataframe = {}
        for k in columns:
            if of_type == 'P':
                column_name = k.split('.')
                if column_name[0] in prop_columns and f.loc[len(f[k])-1, k] == 'SimWell':
                    dataframe[f'{k[:-2]}'] = f[k][:-1]
            else:
                if k in prop_columns:
                    dataframe[f'{k}'] = f[k]
        return dataframe


class PerfilFacies:
    def __init__(self, wellinfo, well_markers, simdir, cell_ref):
        self.name = wellinfo.well.WELL.value
        self.markers_table = well_markers
        self.rowdata = {'depth': wellinfo.curves.DEPT.data, 'facies': wellinfo.curves.fac.data}
        self.sequence = df(self.rowdata)
        self.sequence = Functions.depth_to_thickness(self.sequence)
        self.sequence = self.applyWellMarkers()
        self.root_result = simdir
        self.ref = cell_ref
        self.df_analysis = self.get_well_data()

    def applyWellMarkers(self):
        if type(self.markers_table) != str:
            zones = self.markers_table.columns.to_list()
            idx = self.markers_table.index[self.markers_table['Well'] == self.name]
            top = float(self.markers_table.loc[idx, zones[1]])
            base = float(self.markers_table.loc[idx, zones[2]])

            self.sequence.drop(self.sequence[self.sequence.depth < top].index, inplace=True)
            self.sequence.drop(self.sequence[self.sequence.depth > base].index, inplace=True)
        else:
            pass
        return self.sequence

    def get_well_data(self) -> df:
        file = df(pd.read_excel(self.root_result, sheet_name=f'{self.ref}{self.name}')).set_index('Unnamed: 0')
        list_col = ['depth', 'Thickness', 'facies', 'Thickness.1']
        extract_cols = []
        for col in file.columns:
            if col not in list_col:
                extract_cols.append(col)
            else:
                if file[col].ID == 'GeoWell':
                    extract_cols.append(col)
        file.drop(columns=extract_cols, inplace=True)
        if any(True for col in file.columns if col == 'Thickness.1'):
            file = file.rename(columns={'Thickness.1': 'Thickness'})
        file = file.iloc[:-1]
        return file


class RelatorioPDF(FPDF):
    def __int__(self):
        self.subtitle()

    def prescript(self, modelname):
        self.modelname = modelname

    def header(self):
        logo_folder = os.path.join(os.getcwd(), 'Logo')
        self.image(f'{logo_folder}\\logo_udesc.png', 180, 12, 15)
        self.image(f'{logo_folder}\\logo_forward.png', 12, 13.8, 28)
        self.set_font('times', 'I', size=12)
        self.cell(0, 16, f'Partial Calibration Report - {self.modelname}', ln=True, align='C')
        self.ln(16)
    def subtitle(self, of_type, finaltime, best_variable, OF_value, labels):
        self.set_font('times', '', size=14.8)

        if self.page_no() == 1 and of_type == 'P':
            self.cell(0, 0, 'OF PROOF Calibration Results', align='C')
            self.ln(12)
            self.set_font('times', '', size=12)
            self.cell(0, 0, 'Calibration time: ' + finaltime + 'h', ln=True, align='C')
            self.ln(7)
            self.multi_cell(0, 6.5, 'OF value:  ' + OF_value,
                            align='C')
            self.set_font('times', '', size=10)
            self.ln(7)


            parameters = labels

            best_variable_values = best_variable.split(',')
            col_widths = [22.5] * len(parameters)
            total_width = sum(col_widths)
            page_width = self.w - 2 * self.l_margin
            start_x = (page_width - total_width) / 2 + self.l_margin
            for index, param in enumerate(parameters):
                self.set_font('times', size=8.5)
                self.set_x(start_x + sum(col_widths[:index]))
                truncated_param = (param[:10] + '...') if len(param) > 10 else param  # Truncate titles
                self.cell(col_widths[index], 10, truncated_param, border=1, align='C')
            self.ln()

            for i, value in enumerate(best_variable_values):
                self.set_font('times', size=8.5)
                if i % 2 == 0:
                    self.set_fill_color(173, 216, 230)
                else:
                    self.set_fill_color(255, 255, 255)

                value = value.strip()
                if len(value) > 10:
                    value = value[:10] + '...'

                self.set_x(start_x + sum(col_widths[:i]))
                self.cell(col_widths[i], 10, value, border=1, align='C', fill=True)

            self.ln(12)
        else:
            self.cell(0, 0, 'OF SCOOF Calibration Results', align='C')
            self.ln(12)
            self.set_font('times', '', size=12)
            self.cell(0, 0, 'Calibration time: ' + finaltime + 'h', ln=True, align='C')
            self.ln(7)
            self.multi_cell(0, 6.5, 'OF value:  ' + OF_value,
                            align='C')
            self.set_font('times', '', size=10)
            self.ln(7)

            parameters = labels

            best_variable_values = best_variable.split(',')
            col_widths = [22.5] * len(parameters)
            total_width = sum(col_widths)
            page_width = self.w - 2 * self.l_margin
            start_x = (page_width - total_width) / 2 + self.l_margin
            for index, param in enumerate(parameters):
                self.set_font('times', size=8.5)
                self.set_x(start_x + sum(col_widths[:index]))
                truncated_param = (param[:10] + '...') if len(param) > 10 else param  # Truncate titles
                self.cell(col_widths[index], 10, truncated_param, border=1, align='C')
            self.ln()

            for i, value in enumerate(best_variable_values):
                self.set_font('times', size=8.5)
                if i % 2 == 0:
                    self.set_fill_color(173, 216, 230)
                else:
                    self.set_fill_color(255, 255, 255)

                value = value.strip()
                if len(value) > 10:
                    value = value[:10] + '...'

                self.set_x(start_x + sum(col_widths[:i]))
                self.cell(col_widths[i], 10, value, border=1, align='C', fill=True)

            self.ln(15)
    def footer(self):
        self.set_y(-15)
        self.set_font('times', size=10)
        current_time = datetime.now().strftime('%d-%m-%Y %H:%M:%S')
        self.cell(0, 10, f'Generetad: {current_time} | Page {self.page_no()}', align='C')

    def add_section_title(self, section_name):
        self.ln(5)
        self.set_font('times', 'B', size=12)
        self.cell(0, 10, 'Section: ' + section_name, ln=True, align='L')
        self.ln(5)

    def structure(self, folder_path):

        section_name = os.path.basename(folder_path)

        if section_name == str('Maps'):
            self.add_page()
            self.ln(5)
            self.add_section_title(str('Well maps proportion'))
            self.ln(5)
        if section_name == str('Statistics'):
            self.add_page()
            self.ln(5)
            self.add_section_title('Statistics')
            self.ln(5)
        if section_name == str('P'):
            self.add_page()
            self.ln(5)
            self.add_section_title('PROOF Results Scatter')
            self.ln(5)
        if section_name == str('S'):
            self.add_page()
            self.add_section_title('SCOOF Results Scatter')
            self.ln(5)
        if section_name.endswith('parameters'):
            self.add_section_title('Parameters used in the base model')
            self.ln(5)
        if section_name.endswith('results'):
            self.add_page()
            self.ln(5)
            self.add_section_title('Base model results')
            self.ln(5)

        page_height_limit = 260
        left_margin = 15
        right_margin = 15
        page_width = self.w - right_margin - left_margin

        if not os.path.exists(folder_path):
            self.cell(0, 10, f"Error: The folder '{folder_path}' was not found.", ln=True)
            return

        if not os.path.isdir(folder_path):
            self.cell(0, 10, f"Error: '{folder_path}' is not a valid folder.", ln=True)
            return

        for file in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file)

            try:

                if file.endswith(('jpg', 'png')) and file.startswith('_'):
                    name = os.path.basename(file)

                    if name == str('_@@@structural_parameters.png'):
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Structural Parameters', ln=True, align='C')
                        self.ln(5)

                    if name == str('_eustatic_curve.png'):
                        self.add_page()
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Eustatic Curve', ln=True, align='C')
                        self.ln(8)

                    if name == str('_maps.png'):
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Maps', ln=True, align='C')
                        self.ln(8)

                    if name == str('_Wells.png'):
                        self.add_page()
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Wells', ln=True, align='C')
                        self.ln(5)

                    if name == str('_WellsOBS_Log.png'):
                        self.add_page()
                        self.add_section_title('Well profiles')
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Wells OBS Log', ln=True, align='C')
                        self.ln(5)

                    if name == str('_Wells_Log.png'):
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Wells Log', ln=True, align='C')
                        self.ln(5)

                    if name == str('_Case Graphic - fixed scale in X 0.jpg'):
                        self.add_page()
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Case Graphics', ln=True, align='C')
                        self.ln(5)

                    if self.get_y() + 85 > page_height_limit:
                        self.add_page()
                        self.set_x(left_margin)


                    image_width = page_width - left_margin - right_margin
                    image_height = 95

                    # Center
                    image_x = (page_width - image_width) / 2 + 1.2*left_margin
                    self.image(file_path, image_x, self.get_y(), image_width, image_height)

                    self.ln(image_height + 7)

                elif file.endswith('xlsx') and file.startswith('_'):
                    name = os.path.basename(file)

                    if name == str('_@@sediment_input_parameters.xlsx'):
                        self.add_page()
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Sediment Input Parameters', ln=True, align='C')
                        self.ln(5)

                    if name == str('_sediment_prod_vs_bathymetry.xlsx'):
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Sediment Production Parameters', ln=True, align='C')
                        self.ln(5)

                    if name == str('_transport_coefficients.xlsx'):
                        self.add_page()
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Transport Coefficients Parameters', ln=True, align='C')
                        self.ln(5)

                    if name == str('_Case Statistics.xlsx'):
                        self.add_page()
                        self.ln(5)
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Case Statistics', ln=True, align='C')
                        self.ln(5)

                    if name == str('_Best Simulations.xlsx'):
                        self.set_font('times', size=12)
                        self.cell(0, 0, 'Best Simulations', ln=True, align='C')
                        self.ln(5)


                    df = pd.read_excel(file_path, sheet_name=None)

                    for sheet_name, data in df.items():
                        row_count = len(data)
                        row_height = 10
                        padding = 5
                        total_height = row_count * row_height + (row_count - 1) * padding

                        if self.get_y() + total_height > page_height_limit:
                            self.add_page()
                            self.set_x(left_margin)

                        self.ln(5)

                        # Format Table
                        col_widths = [22.5] * len(data.columns)
                        total_width = sum(col_widths)
                        start_x = (page_width - total_width) / 2 + left_margin  # Centralize

                        # Header
                        for index, header in enumerate(data.columns):
                            self.set_font('times', size=8.5)
                            self.set_x(start_x + sum(col_widths[:index]))
                            truncated_header = (header[:10] + '...') if len(
                                header) > 10 else header  # Truncate titles
                            self.cell(col_widths[index], 10, truncated_header, border=1,
                                      align='C')
                        self.ln()

                        for i, row in data.iterrows():
                            self.set_font('times', size=8.5)
                            if i % 2 == 0:
                                self.set_fill_color(173, 216, 230)  # light blue
                            else:
                                self.set_fill_color(255, 255, 255)  # white

                            for index, item in enumerate(row):
                                self.set_font('times', size=8.5)
                                if index == 0 and len(str(item)) > 10:  # Truncate first column itens
                                    item = (str(item)[:10] + '...')
                                elif index != 0 and isinstance(item, float):  # Rounds items in other columns
                                    item = "{:.5f}".format(item)

                                self.set_x(start_x + sum(col_widths[:index]))
                                self.cell(col_widths[index], 10, str(item), border=1, align='C',
                                          fill=True)

                            self.ln()

                        self.ln(10)

                        if self.get_y() + 10 > page_height_limit:
                            self.add_page()
                            self.set_x(left_margin)

            except Exception as e:

                self.cell(0, 10, f"Error processing {file}: {str(e)}", ln=True)

    def export(self, root_folder, of_type, basemodel, time, parameter, value):

        #workflow_folder = os.path.join(os.getcwd())

        #project_parameters = str(f'{basemodel}'+'_basemodel_parameters')
        #project_results = str(f'{basemodel}'+'_basemodel_results')

        folder_path = [
            f'{root_folder}\\Results\\{of_type}',
            f'{root_folder}\\Results\\{of_type}\\WellPieChart\\Maps',
            f'{root_folder}\\Results\\{of_type}\\Statistics',
        ]

        try:
            self.pdf = RelatorioPDF("P", "mm")
            model_name_str = str(basemodel[5:])
            self.prescript(modelname=model_name_str)
            self.set_auto_page_break(auto=True, margin=15)
            self.add_page()

            if isinstance(value, float):
                OF_value_str = str(value)

            if isinstance(time, float):
                time_str = str(time)

            _, _, labels, _, _ = Statistics.StatisticalAnalysis(root_folder, of_type, param_interval=None).picker()

            best_variable_str = ', '.join(map(str, parameter))
            self.subtitle(of_type, finaltime=time_str, best_variable=best_variable_str, OF_value=OF_value_str, labels=labels[:-1])
            for i in range(len(folder_path)):
                self.structure(folder_path[i])
            time_ = datetime.now().strftime('%d_%m_%Y %Hh%Mm%Ss')

            self.output(f'{root_folder}\\Results\\{of_type}\\Report_{of_type}_'f'{time_}''.pdf')
            print("PDF exported successfully!")

        except Exception as error:
            print("Error when exporting PDF!", {str(error)})