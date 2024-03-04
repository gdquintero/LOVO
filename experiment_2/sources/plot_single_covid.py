import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os
import models

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df_data = pd.read_table(parent+"/data/covid.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_covid.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

# n_train = int(df_data[0].values[0])
# n_test = int(df_data[0].values[1])

# y = np.zeros(n_train)
# y = df_data[0].values[2:len(df_data)-n_test]

# days = [i for i in range(1,n_train+1)]
# x = df_sol[:].values[0]

# t = np.linspace(days[0],days[-1],1000)

# plt.plot(days,y,"o")
# plt.plot(t,models.cubic(x[0],x[1],x[2],t,y[-1],t[-1]))
# plt.show()


# with open(parent+"/output/outliers.txt") as f:
#     lines = f.readlines()
#     xdata = [line.split()[0] for line in lines]

# noutliers = int(xdata[0])
# outliers = np.empty(noutliers,dtype=int)
# cubic_outliers = np.empty((2,noutliers))

# for i in range(noutliers):
#     outliers[i] = int(xdata[i+1])

# for i in range(noutliers):
#     cubic_outliers[0,i] = df_data[0].values[outliers[i]-1]
#     cubic_outliers[1,i] = df_data[1].values[outliers[i]-1]

# t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)

# for i in range(noutliers):
#         point1 = [cubic_outliers[0,i],poly(df_sol.iloc[0].values,cubic_outliers[0,i])]
#         point2 = [cubic_outliers[0,i],cubic_outliers[1,i]]
#         x_values = [point1[0], point2[0]]
#         y_values = [point1[1], point2[1]]
        

# plt.plot(t,poly(df_sol.values[0][:4],t),lw=2)
# plt.plot(cubic_outliers[0],cubic_outliers[1],'ro',mfc='none',ms=10)

# l1 = plt.plot(df_data[0].values[:80],df_data[1].values[:80],"ok")
# plt.setp(l1, 'markersize', 4)

# l2 = plt.plot(df_data[0].values[80:],df_data[1].values[80:],"ok",mfc='none',ms=10,marker="s")
# plt.setp(l2, 'markersize', 4)

# plt.xticks(range(-1, 4, 1))
# plt.yticks(range(-4, 5, 2))
# plt.show()
