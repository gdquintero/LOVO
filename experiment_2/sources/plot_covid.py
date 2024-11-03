import numpy as np
import matplotlib.pyplot as plt
import os
import models
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import pandas as pd

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df_sol_5 = pd.read_table(parent+"/output/plot_covid_5.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_10 = pd.read_table(parent+"/output/plot_covid_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_15 = pd.read_table(parent+"/output/plot_covid_15.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_20 = pd.read_table(parent+"/output/plot_covid_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_25 = pd.read_table(parent+"/output/plot_covid_25.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_data = pd.read_table(parent+"/data/covid_25.txt",delimiter=" ",header=None,skiprows=2,skipinitialspace=True)

t = np.linspace(1,30,30)

plt.plot(df_sol_5[0].values+20,df_sol_5[1].values,label="$\hat m=5$")
plt.plot(df_sol_10[0].values+15,df_sol_10[1].values,label="$\hat m=10$")
plt.plot(df_sol_15[0].values+10,df_sol_15[1].values,label="$\hat m=15$")
plt.plot(df_sol_20[0].values+5,df_sol_20[1].values,label="$\hat m=20$")
plt.plot(df_sol_25[0].values,df_sol_25[1].values,label="$\hat m=25$")

plt.plot(t[:25],df_data[0].values[:25],"ok",alpha=0.7)
plt.plot(t[25:],df_data[0].values[25:],"ok",mfc='none',ms=6,marker="s",alpha=0.7)
# plt.ylim([min(df_data[0].values)-2, max(df_data[0].values)+1])
plt.tick_params(axis='both',direction='in')
plt.legend()
plt.savefig(parent+"/images/image.pdf",bbox_inches = "tight")
plt.show()
