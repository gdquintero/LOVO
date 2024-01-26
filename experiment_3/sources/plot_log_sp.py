import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df_sp = pd.read_table(parent+"/output/log_sp.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

with open(parent+"/output/num_mixed_test.txt") as f:
    lines = f.readlines()
    lim = [line.split()[0] for line in lines]

inf = int(lim[0])
sup = int(lim[1])
n = sup - inf + 1

t = np.linspace(inf,sup,n)

plt.plot(t,np.log(df_sp[0].values),"ko",ls=":")
plt.show()