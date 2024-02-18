import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import models

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df_error = pd.read_table(parent+"/output/error.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

n_train = 30
n_test = 10
samples = 400
ncv = samples // (n_train + n_test)
t = np.linspace(1,ncv,ncv)


plt.ylim(0,1)

plt.rcParams.update({'font.size': 14})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(t,df_error.values[1][:ncv],"-o")
# # plt.plot(t,df_sol_20.values,"-o",label="20")
# # plt.plot(t,df_sol_30.values,"-o",label="30")
# plt.legend()
plt.show()