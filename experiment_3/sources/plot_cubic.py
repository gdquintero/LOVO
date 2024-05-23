import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
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

plt.plot(t,poly(df_sol.values[0][:4],t),lw=1)
plt.plot(cubic_outliers[0],cubic_outliers[1],'ro',mfc='none',ms=6,mew=0.5)
plt.plot(df_data[0].values[:80],df_data[1].values[:80],"ok",ms=2)
plt.plot(df_data[0].values[80:],df_data[1].values[80:],"sk",mfc='none',ms=2)

plt.xticks(range(-1, 4, 1))
plt.yticks(range(-4, 5, 2))
plt.xlabel("$t$")
plt.ylabel("$y(t;x^*)$")
plt.ylim([-4.5,4.5])
plt.tick_params(axis='both',direction='in')
plt.savefig(parent+"/images/cubic.pdf",bbox_inches = "tight")
# plt.show()
