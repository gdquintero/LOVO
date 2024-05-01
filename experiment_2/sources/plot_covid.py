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

df_sol = pd.read_table(parent+"/output/plot_covid.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

plt.plot(df_sol[0].values,df_sol[1].values)
plt.show()
