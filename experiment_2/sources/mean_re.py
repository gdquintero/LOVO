import numpy as np
import pandas as pd
import os

def mount_cumulative_mean():
    cumulative_means[0,0,:] = df_error_5[1].values
    cumulative_means[1,0,:] = df_error_10[1].values
    cumulative_means[2,0,:] = df_error_15[1].values
    cumulative_means[3,0,:] = df_error_20[1].values
    cumulative_means[4,0,:] = df_error_25[1].values
    cumulative_means[5,0,:] = df_error_30[1].values

    for i in range(1,10):
        for j in range(10):
            cumulative_means[0,i,j] = np.mean(df_error_5.values[i][1:j+2])
            cumulative_means[1,i,j] = np.mean(df_error_10.values[i][1:j+2])
            cumulative_means[2,i,j] = np.mean(df_error_15.values[i][1:j+2])
            cumulative_means[3,i,j] = np.mean(df_error_20.values[i][1:j+2])
            cumulative_means[4,i,j] = np.mean(df_error_25.values[i][1:j+2])
            cumulative_means[5,i,j] = np.mean(df_error_30.values[i][1:j+2])

    with open(parent+"/output/cumulative_mean_re_5.txt","w") as f:
        for i in range(10):
            f.write("%f %f %f %f %f %f %f %f %f %f\n" % (
                cumulative_means[0,i,0],cumulative_means[0,i,1],cumulative_means[0,i,2],cumulative_means[0,i,3],cumulative_means[0,i,4],\
                cumulative_means[0,i,5],cumulative_means[0,i,6],cumulative_means[0,i,7],cumulative_means[0,i,8],cumulative_means[0,i,9]))
            
        with open(parent+"/output/cumulative_mean_re_10.txt","w") as f:
            for i in range(10):
                f.write("%f %f %f %f %f %f %f %f %f %f\n" % (
                    cumulative_means[1,i,0],cumulative_means[1,i,1],cumulative_means[1,i,2],cumulative_means[1,i,3],cumulative_means[1,i,4],\
                    cumulative_means[1,i,5],cumulative_means[1,i,6],cumulative_means[1,i,7],cumulative_means[1,i,8],cumulative_means[1,i,9]))
                
        with open(parent+"/output/cumulative_mean_re_15.txt","w") as f:
            for i in range(10):
                f.write("%f %f %f %f %f %f %f %f %f %f\n" % (
                    cumulative_means[2,i,0],cumulative_means[2,i,1],cumulative_means[2,i,2],cumulative_means[2,i,3],cumulative_means[2,i,4],\
                    cumulative_means[2,i,5],cumulative_means[2,i,6],cumulative_means[2,i,7],cumulative_means[2,i,8],cumulative_means[2,i,9]))
                
        with open(parent+"/output/cumulative_mean_re_20.txt","w") as f:
            for i in range(10):
                f.write("%f %f %f %f %f %f %f %f %f %f\n" % (
                    cumulative_means[3,i,0],cumulative_means[3,i,1],cumulative_means[3,i,2],cumulative_means[3,i,3],cumulative_means[3,i,4],\
                    cumulative_means[3,i,5],cumulative_means[3,i,6],cumulative_means[3,i,7],cumulative_means[3,i,8],cumulative_means[3,i,9]))
                
        with open(parent+"/output/cumulative_mean_re_25.txt","w") as f:
            for i in range(10):
                f.write("%f %f %f %f %f %f %f %f %f %f\n" % (
                    cumulative_means[4,i,0],cumulative_means[4,i,1],cumulative_means[4,i,2],cumulative_means[4,i,3],cumulative_means[4,i,4],\
                    cumulative_means[4,i,5],cumulative_means[4,i,6],cumulative_means[4,i,7],cumulative_means[4,i,8],cumulative_means[4,i,9]))
                
        with open(parent+"/output/cumulative_mean_re_30.txt","w") as f:
            for i in range(10):
                f.write("%f %f %f %f %f %f %f %f %f %f\n" % (
                    cumulative_means[5,i,0],cumulative_means[5,i,1],cumulative_means[5,i,2],cumulative_means[5,i,3],cumulative_means[5,i,4],\
                    cumulative_means[5,i,5],cumulative_means[5,i,6],cumulative_means[5,i,7],cumulative_means[5,i,8],cumulative_means[5,i,9]))
                
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

mount_cumulative_mean()
mount_mean_per_day()
