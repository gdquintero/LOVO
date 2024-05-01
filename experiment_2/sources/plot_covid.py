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

df_data = pd.read_table(parent+"/data/covid.txt",delimiter=" ",header=None,skiprows=2,skipinitialspace=True)
df_sols = pd.read_table(parent+"/output/solutions_5_25.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)



plt.show()
