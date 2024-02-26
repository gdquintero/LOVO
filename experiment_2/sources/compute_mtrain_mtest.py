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

successfull_threshold = np.zeros((6,10),dtype=int)

for i in range(10):
    for j in range(10):
        if df_error_5[i+1].values[j] <= threshold:
            successfull_threshold[0,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_10[i+1].values[j] <= threshold:
            successfull_threshold[1,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_15[i+1].values[j] <= threshold:
            successfull_threshold[2,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_20[i+1].values[j] <= threshold:
            successfull_threshold[3,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_25[i+1].values[j] <= threshold:
            successfull_threshold[4,i] += 1

for i in range(10):
    for j in range(10):
        if df_error_30[i+1].values[j] <= threshold:
            successfull_threshold[5,i] += 1
            

with open(parent+"/output/successfull_threshold.txt","w") as f:
        for i in range(6):
            f.write("%i %i %i %i %i %i %i %i %i %i\n" % (
                    successfull_threshold[i,0],successfull_threshold[i,1],successfull_threshold[i,2],successfull_threshold[i,3],successfull_threshold[i,4],\
                    successfull_threshold[i,5],successfull_threshold[i,6],successfull_threshold[i,7],successfull_threshold[i,8],successfull_threshold[i,9]))