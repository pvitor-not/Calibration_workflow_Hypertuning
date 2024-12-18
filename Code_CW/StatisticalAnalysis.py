import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import seaborn as sns
from scipy.stats import hmean, gmean
from Code_CW import Functions
from Code_CW.Objects import Simulation_parameters


class StatisticalAnalysis:
    def __init__(self, root_folder, of_type, param_interval=None):
        self.root = root_folder
        self.of_type = of_type.upper()
        self.picker()
        self.statistics()
        self.plot(param_interval)
        # self.dionisosparam()

    def picker(self):
        try:
            os.makedirs(f'{self.root}\\Results\\{self.of_type}\\Statistics', exist_ok=True)
            file_path = f'{self.root}\\Results\\{self.of_type}\\Uncertain_parameters_and_OF_values.xlsx'
            xlsx_best_path = f'{self.root}\\Results\\{self.of_type}\\Statistics\\_Best Simulations.xlsx'

            simdata = pd.read_excel(file_path)
            simdata = simdata.loc[:, ~simdata.columns.isin(['Ambiguity', 'Unnamed: 0', 'Simulation', 'OutputPath'])]
            lenght = len(simdata.index)

            percent = 15

            lenght2 = round(percent * 1 / 100 * lenght)
            simdata = simdata.sort_values(by='OF Value')
            labels = list(simdata.columns)

            simdata_best = simdata.iloc[:lenght2, :]  # Remove the lines from the datalen2 to the last one
            simdata_best.to_excel(xlsx_best_path)
            pd.options.display.float_format = '{:.5f}'.format  # Verificar aqui

            return simdata, simdata_best, labels, lenght, lenght2

        except Exception as e:
            print("Error in statistics function!")
            print(f"An error occurred in picker: {e}")

    def statistics(self):
        try:
            xlsx_stats_path = f'{self.root}\\Results\\{self.of_type}\\Statistics\\_Case Statistics.xlsx'

            simdata, simdata_best, labels, datalen, datalen2 = self.picker()

            maxi = simdata_best.max()
            mini = simdata_best.min()
            mean = simdata_best.mean()
            median = simdata_best.median()
            geo_mean = gmean(simdata_best)
            geo_mean = pd.DataFrame({"Column": np.array(geo_mean)})
            geo_mean.index = labels
            harmonic_mean = hmean(simdata_best)
            harmonic_mean = pd.DataFrame({"Column": np.array(harmonic_mean)})
            harmonic_mean.index = labels
            sample_std = simdata_best.std(ddof=1)

            index = ['Máximo', 'Mínimo', 'Média', 'Mediana', 'Média Geométrica', 'Média Harmônica',
                     'Desvio Padrão Amostral']

            df_stats = [maxi, mini, mean, median, geo_mean['Column'], harmonic_mean['Column'], sample_std]

            df_stats = pd.DataFrame(df_stats, index=index, columns=labels)

            df_stats.to_excel(xlsx_stats_path)

        except Exception as e:
            print("Error in statistics function!")
            print(f"An error occurred: {e}")
        else:
            print("Multipliers statistics successful done!")
            return simdata, simdata_best, labels, datalen, datalen2

    def plot(self, param_interval):
        try:
            if param_interval is None:
                param_interval = [0.5, 1.5]  # Global
            bars_qnt = round((np.array(param_interval).max() - np.array(param_interval).min()) / 0.1)
            simdata, simdata_best, labels, datalen, datalen2 = self.statistics()
            config = ['fixed scale in X', 'non-fixed scale in X']

            for k in range(len(config)):
                remain_var = len(labels) - 1  # excluding the OF Values
                it = 0
                while remain_var >= 1:

                    # Set graphics size
                    if remain_var >= 4:
                        j = 4
                    elif remain_var == 3:
                        j = 3
                    elif remain_var == 2:
                        j = 2
                    else:
                        j = 1

                    sns.set(font_scale=0.36)
                    fig, axs = plt.subplots(nrows=2, ncols=j)

                    for i in range(j):

                        if config[k] == 'fixed scale in X':
                            axs[0, i].set_xlim(param_interval)
                            axs[1, i].set_xlim(param_interval)
                        axs[0, i].set_ylim([0, 100])
                        axs[1, i].set_ylim([0, 100])
                        axs[0, i].set_ylabel(" ")
                        axs[0, i].set_xlabel(" ")
                        axs[1, i].set_ylabel(" ")
                        axs[0, i].yaxis.set_major_formatter(PercentFormatter(100))
                        axs[1, i].yaxis.set_major_formatter(PercentFormatter(100))
                        if i != 0:
                            axs[0, i].axes.yaxis.set_ticklabels([])
                            axs[1, i].axes.yaxis.set_ticklabels([])

                        if i == 0 and labels[4 * it] not in "OF Values":

                            sns.histplot(x=labels[4 * it], data=simdata, bins=bars_qnt, stat="percent", kde=False,
                                         ax=axs[0, 0], binrange=param_interval, legend=False, color=[0.1, 0, 0.3, 0.8])
                            sns.histplot(x=labels[4 * it], data=simdata_best, bins=bars_qnt, stat="percent", kde=False,
                                         ax=axs[1, 0], binrange=param_interval, legend=False, color=[0.1, 0, 0.3, 0.8])
                        elif i == 1 and labels[4 * it + 1] not in "OF Values":

                            sns.histplot(x=labels[4 * it + 1], data=simdata, bins=bars_qnt, stat="percent", kde=False,
                                         ax=axs[0, 1], binrange=param_interval, legend=False,
                                         color=[0.4, 0.5, 0.2, 0.8])
                            sns.histplot(x=labels[4 * it + 1], data=simdata_best, bins=bars_qnt, stat="percent",
                                         kde=False,
                                         ax=axs[1, 1], binrange=param_interval, legend=False,
                                         color=[0.4, 0.5, 0.2, 0.8])
                        elif i == 2 and labels[4 * it + 2] not in "OF Values":

                            sns.histplot(x=labels[4 * it + 2], data=simdata, bins=bars_qnt, stat="percent", kde=False,
                                         ax=axs[0, 2], binrange=param_interval, legend=False,
                                         color=[0.1, 0.3, 0.5, 0.8])
                            sns.histplot(x=labels[4 * it + 2], data=simdata_best, bins=bars_qnt, stat="percent",
                                         kde=False,
                                         ax=axs[1, 2], binrange=param_interval, legend=False,
                                         color=[0.1, 0.3, 0.5, 0.8])
                        elif i == 3 and labels[4 * it + 3] not in "OF Values":
                            sns.histplot(x=labels[4 * it + 3], data=simdata, bins=bars_qnt, stat="percent", kde=False,
                                         ax=axs[0, 3], binrange=param_interval, legend=False,
                                         color=[0.3, 0.2, 0.1, 0.8])
                            sns.histplot(x=labels[4 * it + 3], data=simdata_best, bins=bars_qnt, stat="percent",
                                         kde=False,
                                         ax=axs[1, 3], binrange=param_interval, legend=False,
                                         color=[0.3, 0.2, 0.1, 0.8])

                    title = "Case analyses" + "\n" + "Number of simulations: " + str(
                        datalen) + "     " + "Best simulations: " + str(
                        datalen2) + "     " + "Bars length: " + str(
                        (np.array(param_interval).max() - np.array(param_interval).min()) / bars_qnt)
                    fig.suptitle(title)

                    plt.subplots_adjust(left=0.043, bottom=0.235, right=0.99, top=0.926, wspace=0.180, hspace=0.130)

                    png_path = f'{self.root}\\Results\\{self.of_type}\\Statistics\\_Case Graphic - {str(config[k])} {str(it)}.jpg'

                    fig.savefig(png_path, format='jpg', dpi=800)

                    remain_var = remain_var - 4
                    it += 1

        except Exception as e:
            print("Error in plotter function!")
            print(f"An error occurred: {e}")

        else:
            print("Plotting done!")
