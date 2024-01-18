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
df_sol = pd.read_table(parent+"/output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

with open(parent+"/output/outliers.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

noutliers = int(xdata[0])
outliers = np.empty(noutliers,dtype=int)
cubic_outliers = np.empty((2,noutliers))

for i in range(noutliers):
    outliers[i] = int(xdata[i+1])

for i in range(noutliers):
    cubic_outliers[0,i] = df_data[0].values[outliers[i]-1]
    cubic_outliers[1,i] = df_data[1].values[outliers[i]-1]

t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)

# plt.plot(df_data[0].values,df_data[1].values,"ok")
plt.plot(t,poly(df_sol.values[0][:4],t))

for i in range(noutliers):
        point1 = [cubic_outliers[0,i],poly(df_sol.iloc[0].values,cubic_outliers[0,i])]
        point2 = [cubic_outliers[0,i],cubic_outliers[1,i]]
        x_values = [point1[0], point2[0]]
        y_values = [point1[1], point2[1]]
        # plt.plot(x_values, y_values, 'k', linestyle="--")
        
l = plt.plot(df_data[0].values,df_data[1].values,"ko")
plt.setp(l, 'markersize', 6)
    
plt.show()
