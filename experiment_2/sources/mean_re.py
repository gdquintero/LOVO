import numpy as np
import pandas as pd
import os
                
def mount_mean_per_day():
    for i in range(10):
        mean_per_day[0,i] = np.mean(df_error_5[i+1].values)
        mean_per_day[1,i] = np.mean(df_error_10[i+1].values)
        mean_per_day[2,i] = np.mean(df_error_15[i+1].values)
        mean_per_day[3,i] = np.mean(df_error_20[i+1].values)
        mean_per_day[4,i] = np.mean(df_error_25[i+1].values)
        mean_per_day[5,i] = np.mean(df_error_30[i+1].values)

    with open(parent+"/output/mean_per_day_re.txt","w") as f:
        for i in range(6):
            f.write("%f %f %f %f %f %f %f %f %f %f\n" % (
                    mean_per_day[i,0],mean_per_day[i,1],mean_per_day[i,2],mean_per_day[i,3],mean_per_day[i,4],\
                    mean_per_day[i,5],mean_per_day[i,6],mean_per_day[i,7],mean_per_day[i,8],mean_per_day[i,9]))


cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df_error_5 = pd.read_table(parent+"/output/latex_5.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_10 = pd.read_table(parent+"/output/latex_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_15 = pd.read_table(parent+"/output/latex_15.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_20 = pd.read_table(parent+"/output/latex_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_25 = pd.read_table(parent+"/output/latex_25.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_30 = pd.read_table(parent+"/output/latex_30.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

cumulative_means = np.zeros((6,10,10))    
mean_per_day = np.zeros((6,10))

mount_mean_per_day()
