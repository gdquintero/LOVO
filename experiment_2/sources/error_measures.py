import numpy as np
import pandas as pd
import os
import models
import statistics as st
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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

for i in range(5):
    print("%.4f" % st.mean(errors[:,i]))

successes = np.zeros(5,dtype=int)

for i in range(1000):
    for j in range(5):
        if errors[i,j] <= 0.1:
            successes[j] += 1


print("\nSuccesses: ")
print(successes)

fig, ax = plt.subplots()
mtrains = df_optimal_ntrains[0].values[:]
labels = ["%s"%i for i in range(5,31,5)]
unique, counts = np.unique(mtrains, return_counts=True)
rects = ax.bar(unique,counts, 3)
ax.set_xticks(unique)
ax.set_xticklabels(labels)

print("\nFrecuencia: ")
print(counts/10)

plt.xlabel('Days for training',fontsize=16)
plt.ylabel('Frequency',fontsize=16)
plt.show()