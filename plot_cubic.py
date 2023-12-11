import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def poly(x,t):
    return x[0] + x[1] * t + x[2] * (t**2) + x[3] * (t**3)

df_data = pd.read_table("data/cubic.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)
df_sol = pd.read_table("output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)

plt.plot(df_data[0].values,df_data[1].values,"o")

# for i in range(len(df_sol)):
#     plt.plot(t,poly(df_sol.values[i][:4],t))


plt.plot(t,poly(df_sol.values[0][:4],t))
plt.plot(t,poly(df_sol.values[5][:4],t))

plt.show()
