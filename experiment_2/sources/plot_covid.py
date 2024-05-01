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

x = df_sols[:].values[:]
y = df_data[0].values[:]
days = [i for i in range(1,31)]
t = np.linspace(days[0],days[-1],1000)

plt.plot(t,models.cubic(x[4,0],x[4,1],x[4,2],t,y[len(y)-6],days[len(y)-6]),lw=2)

l1 = plt.plot(days[:25],y[:25],"ok")
plt.setp(l1, 'markersize')

l2 = plt.plot(days[25:],y[25:],"ok",mfc='none',ms=6,marker="s")
plt.setp(l2, 'markersize')

plt.show()
