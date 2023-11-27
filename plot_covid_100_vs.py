import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import models

df_sol_cubic = pd.read_table("output/mean_percentage_error_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_logistic = pd.read_table("output/mean_percentage_error_logistic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

t = np.linspace(1,100,100)

plt.ylim(0,100)

plt.rcParams.update({'font.size': 14})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(t,df_sol_cubic.values,"-o",label="cubic")
plt.plot(t,df_sol_logistic.values,"-o",label="logistic")
plt.legend()
plt.show()