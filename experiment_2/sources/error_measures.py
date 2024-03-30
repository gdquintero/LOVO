import numpy as np
import pandas as pd
import os
import models
import statistics as st

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

def rmsd (n,o,p):

    res = 0

    for i in range(n):
        res = res + (o[i]-p[i])**2

    return np.sqrt(res / n)

def relative_error(o,p):
    return abs(p - o) / abs(o)

df_data = pd.read_table(parent+"/data/covid_mixed.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)
df_sols = pd.read_table(parent+"/output/solutions_covid_mixed.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_optimal_ntrains = pd.read_table(parent+"/output/optimal_ntrains.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

errors = np.empty((1000,5))

for i in range(1000):
    x = df_sols[:].values[i]
    tm = df_optimal_ntrains[0].values[i]
    ym = df_data[0].values[i+34]

    for j in range(5):
        p = models.cubic(x[0],x[1],x[2],tm+j+1,ym,tm)
        o = df_data[0].values[i+j+35]
        errors[i,j] = relative_error(o,p)


print("\nAverages: ")
print(st.mean(errors[:,0]))
print(st.mean(errors[:,1]))
print(st.mean(errors[:,2]))
print(st.mean(errors[:,3]))
print(st.mean(errors[:,4]))

successes = np.zeros(5,dtype=int)

for i in range(1000):
    for j in range(5):
        if errors[i,j] <= 0.1:
            successes[j] += 1

print(successes)
