import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import models
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))


def plot_models(opt=None):
    plt.rcParams.update({'font.size': 14})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylim(min(min(y[n_train - sup:]),min(y_later))-0.5,max(max(y[n_train - sup:]),max(y_later))+0.5)


    plt.plot(t,models.cubic(*df_solution_cubic.values[0][:],t,y[-1],days[-1]))
    plt.plot(outliers_covid[0,:],outliers_covid[1,:],'ro',mfc='none',ms=10)
    # plt.savefig("images/single_covid_cubic.pdf",bbox_inches = "tight")

    plt.plot(np.linspace(1,sup,sup),y[len(y)-sup:],"ko")
    # plt.plot(later_days,y_later,"^",color="grey")
    plt.plot(days_later,y_later,"k^",mfc='none')
    plt.show()
    plt.close()

# Reading the necessary data
df_solution_cubic       = pd.read_table(parent+"/output/solutions_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_train_set            = pd.read_table(parent+"/data/covid_train.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_test_set             = pd.read_table(parent+"/data/covid_test.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
n_train                 = int(df_train_set.values[0][0])
n_test                  = int(df_test_set.values[0][0])

with open(parent+"/output/inf_sup_covid.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

sup = int(xdata[1])

# Arrays allocation
y               = np.zeros(n_train)
y_later         = np.zeros(n_test)
days            = np.linspace(1,sup,sup)
days_later      = np.linspace(1 + sup,sup + n_test,n_test)
t               = np.linspace(1,sup + n_test,1000)

# Observation in the previous days considered
for i in range(n_train):
    y[i] = df_train_set.values[i+1][0]

# Observation in the following days
for i in range(n_test):
    y_later[i] = df_test_set.values[i+1][0]

with open(parent+"/output/outliers_covid_cubic.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

noutliers_covid     = int(xdata[0])
outliers_covid      = np.empty((2,noutliers_covid))
ind_outliers_covid  = np.empty(noutliers_covid,dtype=int)

for i in range(noutliers_covid):
    ind_outliers_covid[i] = int(xdata[i+1])

for i in range(noutliers_covid):
    outliers_covid[0,i] = days[ind_outliers_covid[i]-1]
    outliers_covid[1,i] = y[n_train - sup + ind_outliers_covid[i]-1]


# with open("output/outliers_covid_der_logistic.txt") as f:
#     lines = f.readlines()
#     xdata = [line.split()[0] for line in lines]

# for i in range(noutliers):
#     ind_outliers[i] = int(xdata[i+1])

# for i in range(noutliers):
#     outliers[2,0,i] = days[ind_outliers[i]-1]
#     outliers[2,1,i] = y[n_train - sup + ind_outliers[i]-1]

plot_models(1)
# plot_models(2)
# plot_models()
