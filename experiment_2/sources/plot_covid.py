import numpy as np
import matplotlib.pyplot as plt
import os
import models
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import pandas as pd

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df_sol_5 = pd.read_table(parent+"/output/plot_covid_5.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_10 = pd.read_table(parent+"/output/plot_covid_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_15 = pd.read_table(parent+"/output/plot_covid_15.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_20 = pd.read_table(parent+"/output/plot_covid_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_25 = pd.read_table(parent+"/output/plot_covid_25.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_data = pd.read_table(parent+"/data/covid_25.txt",delimiter=" ",header=None,skiprows=2,skipinitialspace=True)

t = np.linspace(1,30,30)

plt.plot(df_sol_5[0].values+20,df_sol_5[1].values)
plt.plot(df_sol_10[0].values+15,df_sol_10[1].values)
plt.plot(df_sol_15[0].values+10,df_sol_15[1].values)
plt.plot(df_sol_20[0].values+5,df_sol_20[1].values)
plt.plot(df_sol_25[0].values,df_sol_25[1].values)

l1 = plt.plot(t[:25],df_data[0].values[:25],"ok")
plt.setp(l1, 'markersize')

l2 = plt.plot(t[25:],df_data[0].values[25:],"ok",mfc='none',ms=6,marker="s")
plt.setp(l2, 'markersize')
plt.savefig(parent+"/images/image.pdf",bbox_inches = "tight")
# plt.show()
