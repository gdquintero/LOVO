import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def cubic(x1,x2,x3,t,ym,tm):
    return ym + x1 * (t - tm) + x2 * (t - tm)**2 + x3 * (t - tm)**3

df_solution = pd.read_table("output/solutions_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_train_set = pd.read_table("output/covid_train.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
n_train = int(df_train_set.values[0][0])

x   = np.zeros(3)
y   = np.zeros(n_train)
day = np.linspace(1,20,20)
t   = np.linspace(1,20,1000)

x[0] = df_solution.values[0][0]
x[1] = df_solution.values[0][1]
x[2] = df_solution.values[0][2]

for i in range(n_train):
    y[i] = df_train_set.values[i+1][0]

with open("output/outliers_covid_cubic.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

noutliers = int(xdata[0])

ind_outliers = np.empty(noutliers,dtype=int)

for i in range(noutliers):
    ind_outliers[i] = int(xdata[i+1])

outliers = np.empty((2,noutliers))

for i in range(noutliers):
    outliers[0,i] = day[ind_outliers[i]-1]
    outliers[1,i] = y[n_train - 20 + ind_outliers[i]-1]


# plt.plot(day,y[n_train-20:],"ko")
plt.plot(t,cubic(*x,t,y[-1],t[-1]))
plt.plot(outliers[0],outliers[1],'ro',mfc='none',ms=10)

for i in range(noutliers):
    point1 = [outliers[0,i],cubic(*x,outliers[0,i],y[-1],t[-1])]
    point2 = [outliers[0,i],outliers[1,i]]
    x_values = [point1[0], point2[0]]
    y_values = [point1[1], point2[1]]
    plt.plot(x_values, y_values, 'k', linestyle="--")



l = plt.plot(day,y[n_train-20:],"ko")
plt.setp(l, 'markersize', 6)

plt.savefig("cubic.pdf",bbox_inches = "tight")
plt.show()
