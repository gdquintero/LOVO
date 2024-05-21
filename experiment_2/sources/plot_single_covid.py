import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import models

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df_data = pd.read_table(parent+"/data/covid.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_covid.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

n_train = int(df_data[0].values[0])
n_test = int(df_data[0].values[1])

y = df_data[0].values[2:len(df_data)]

days = [i for i in range(1,n_train+n_test+1)]
x = df_sol[:].values[0]

t = np.linspace(days[0],days[-1],1000)

with open(parent+"/output/outliers.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

noutliers = int(xdata[0])
outliers = np.empty(noutliers,dtype=int)
cubic_outliers = np.empty((2,noutliers))

for i in range(noutliers):
    outliers[i] = int(xdata[i+1])

for i in range(noutliers):
    cubic_outliers[0,i] = days[outliers[i]-1]
    cubic_outliers[1,i] = df_data.values[outliers[i]+1]

plt.plot(t,models.cubic(x[0],x[1],x[2],t,y[len(y)-n_test-1],days[len(y)-n_test-1]),lw=1)
plt.plot(cubic_outliers[0],cubic_outliers[1],'ro',mfc='none',ms=6)

plt.plot(days[:n_train],y[:n_train],"ok",ms=3)
plt.plot(days[n_train:],y[n_train:],"sk",mfc='none',ms=3)

plt.xticks(np.arange(0,31,5))
plt.yticks(np.arange(0.1,0.55,0.1))
plt.xlim(-1,31)
plt.ylim(0.05,0.5)
plt.tick_params(axis='both',direction='in')
plt.savefig(parent+"/images/image.pdf",bbox_inches = "tight")
# plt.show()
