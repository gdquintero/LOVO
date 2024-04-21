import numpy as np
import pandas as pd
import os
import models
import statistics as st
import matplotlib.pyplot as plt

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))
population = 41128772

def rmsd (n,o,p):

    res = 0

    for i in range(n):
        res = res + (o[i]-p[i])**2

    return np.sqrt(res / n)

def relative_error(o,p):
    return abs(p - o) / abs(o)

def absolute_error(o,p):
    return abs(o-p)

df_data = pd.read_table(parent+"/data/covid.txt",delimiter=" ",header=None,skiprows=2,skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_covid.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

tm = 27
ym = df_data[0].values[26]
x = df_sol[:].values[0]

errors = np.empty(3)
pred = np.empty(3)
obs = np.empty(3)

for i in range(3):
    pred[i] = models.cubic(x[0],x[1],x[2],tm+i+1,ym,tm)
    obs[i] = df_data[0].values[i]
    
errors[:] = absolute_error(obs,pred)

print("Obs: ",obs)

print("\nPredictions: ")
for i in range(3):
    print("%.3f" % pred[i])

print("\nAbsolute errors: ")
for i in range(3):
    print("%.3f" % errors[i])

print("\nMean AE: ")
print("%.3f" % st.mean(errors))

# print("\nRMSD: ")
# print("%.4f" % rmsd(3,obs,pred))

