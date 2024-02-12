import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import models

df_sol_10 = pd.read_table("output/mean_percentage_error_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_15 = pd.read_table("output/mean_percentage_error_15.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_20 = pd.read_table("output/mean_percentage_error_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_25 = pd.read_table("output/mean_percentage_error_25.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol_30 = pd.read_table("output/mean_percentage_error_30.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

t = np.linspace(1,100,100)

plt.ylim(0,100)

plt.rcParams.update({'font.size': 14})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(t,df_sol_10.values,"-o",label="10")
plt.plot(t,df_sol_15.values,"-o",label="15")
# plt.plot(t,df_sol_20.values,"-o",label="20")
# plt.plot(t,df_sol_25.values,"-o",label="25")
# plt.plot(t,df_sol_30.values,"-o",label="30")
plt.legend()
plt.show()