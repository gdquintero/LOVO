import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def cubic(x1,x2,x3,t,ym,tm):
    return ym + x1 * (t - tm) + x2 * (t - tm)**2 + x3 * (t - tm)**3

df_solution = pd.read_table("output/solution_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_train_set = pd.read_table("output/covid_train.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
n_train = int(df_train_set.values[0][0])

x = np.zeros(3)
y = np.zeros(n_train)
day = np.linspace(1,10,10)
t = np.linspace(1,10,1000)

x[0] = df_solution.values[0][0]
x[1] = df_solution.values[0][1]
x[2] = df_solution.values[0][2]

for i in range(n_train):
    y[i] = df_train_set.values[i+1][0]


plt.plot(day,y[n_train-10:],"ko")
plt.plot(t,cubic(*x,t,y[-1],10))
plt.savefig("cubic.pdf",bbox_inches = "tight")
plt.close()