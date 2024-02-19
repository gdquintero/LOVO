import numpy as np
import pandas as pd
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df_error_10 = pd.read_table(parent+"/output/latex_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_20 = pd.read_table(parent+"/output/latex_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_error_30 = pd.read_table(parent+"/output/latex_30.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

means = np.zeros((3,10,10))
means[0,0,:] = df_error_10[0].values
means[1,0,:] = df_error_20[0].values
means[2,0,:] = df_error_30[0].values

for i in range(10):
    for j in range(10):
        means[0,i,j] = np.mean(df_error_10.values[i][1:j+2])
        means[1,i,j] = np.mean(df_error_20.values[i][1:j+2])
        means[2,i,j] = np.mean(df_error_30.values[i][1:j+2])


# with open(parent+"/output/mean_re_10.txt","w") as f:
#     for i in range(10):
#         f.write("%f %f %f %f %f %f %f %f %f %f\n" % (means[0][0],means[0][1],means[0][2],means[0][3],means[0][4],\
#                                                      means[0][5],means[0][6],means[0][7],means[0][8],means[0][9]))
    
# with open(parent+"/output/mean_re_20.txt","w") as f:
#     f.write("%f %f %f %f %f %f %f %f %f %f\n" % (means[1][0],means[1][1],means[1][2],means[1][3],means[1][4],\
#                                                  means[1][5],means[1][6],means[1][7],means[1][8],means[1][9]))
    
# with open(parent+"/output/mean_re_30.txt","w") as f:
#     f.write("%f %f %f %f %f %f %f %f %f %f\n" % (means[2][0],means[2][1],means[2][2],means[2][3],means[2][4],\
#                                                  means[2][5],means[2][6],means[2][7],means[2][8],means[2][9]))
        

