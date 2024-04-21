import numpy as np
import pandas as pd
import os
import models
import statistics as st
import matplotlib.pyplot as plt

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

def absolute_error(o,p):
    return abs(o-p)

def poly(x,t):
    return x[0] + x[1] * t + x[2] * (t**2) + x[3] * (t**3)

df_data = pd.read_table(parent+"/data/cubic.txt",delimiter=" ",header=None,skiprows=2,skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

print(df_data)

# tm = 27
# ym = df_data[0].values[26]
# x = df_sol[:].values[0]

# errors = np.empty(3)
# pred = np.empty(3)
# obs = np.empty(3)

# # Absolute error in training data

# with open(parent+"/output/outliers.txt") as f:
#     lines = f.readlines()
#     xdata = [line.split()[0] for line in lines]

# noutliers = int(xdata[0])
# outliers = np.empty(noutliers,dtype=int)
# cubic_outliers = np.empty((2,noutliers))

# for i in range(noutliers):
#     outliers[i] = int(xdata[i+1]) - 1

# maeo = 0

# for i in range(tm - noutliers):
#     if i not in outliers:
#         maeo += abs(df_data[0].values[i] - models.cubic(x[0],x[1],x[2],i,ym,tm))
    
# maeo /= (tm - noutliers)

# print("\nMAEO: ")
# print("%.3f" % maeo)
