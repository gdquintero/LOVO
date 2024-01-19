import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def poly(x,t):
    return x[0] + x[1] * t + x[2] * (t**2) + x[3] * (t**3)

df_data = pd.read_table(parent+"/data/cubic.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)

t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)
l = plt.plot(df_data[0].values,df_data[1].values,"ok")
plt.setp(l, 'markersize', 4)
plt.show()