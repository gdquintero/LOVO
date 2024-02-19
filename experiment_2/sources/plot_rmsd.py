import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import models

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df_error_10 = pd.read_table(parent+"/output/error_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_20 = pd.read_table(parent+"/output/error_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_30 = pd.read_table(parent+"/output/error_30.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

ncv = 10
t = np.linspace(1,ncv,ncv)
rmsd = np.zeros((3,ncv))

rmsd[0,:] = df_error_10[0].values
rmsd[1,:] = df_error_20[0].values
rmsd[2,:] = df_error_30[0].values

# plt.ylim(0,1.1)

plt.rcParams.update({'font.size': 14})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for i in range(3):
    plt.plot(t,rmsd[i,:],"--o",label="$n:$ " + str((i+1) * 10))

plt.legend()
plt.show()