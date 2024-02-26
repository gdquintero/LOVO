import numpy as np
import pandas as pd
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df_error_5 = pd.read_table(parent+"/output/latex_5.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_10 = pd.read_table(parent+"/output/latex_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_15 = pd.read_table(parent+"/output/latex_15.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_20 = pd.read_table(parent+"/output/latex_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_25 = pd.read_table(parent+"/output/latex_25.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_30 = pd.read_table(parent+"/output/latex_30.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

threshold = 0.2

successfull_treshold = np.zeros((6,10),dtype=int)

for i in range(10):
    for j in range(10):
        if df_error_5[i+1].values[j] <= threshold:
            successfull_treshold[0,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_10[i+1].values[j] <= threshold:
            successfull_treshold[1,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_15[i+1].values[j] <= threshold:
            successfull_treshold[2,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_20[i+1].values[j] <= threshold:
            successfull_treshold[3,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_25[i+1].values[j] <= threshold:
            successfull_treshold[4,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_30[i+1].values[j] <= threshold:
            successfull_treshold[5,i] += 1
            

