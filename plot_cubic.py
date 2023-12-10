import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def poly(x,t):
    return x[0] + x[1] * (t - 3.) + x[2] * ((t - 3.)**2) + x[3] * ((t - 3.)**3)

df_data = pd.read_table("data/cubic.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)
df_sol = pd.read_table("output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

t = np.linspace(0,6,1000)

plt.plot(df_data[0].values,df_data[1].values,"o")
plt.plot(t,poly(*df_sol.values,t))
plt.show()