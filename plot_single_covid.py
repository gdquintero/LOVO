import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import models


def plot_models(opt=None):
    plt.plot(days,y[n_train - previous_days:],"ko")
    plt.plot(days_later,y_later,"o",color="grey")
    plt.ylim(min(min(y[n_train - previous_days:]),min(y_later))-0.5,max(max(y[n_train - previous_days:]),max(y_later))+0.5)

    if opt == 1:
        plt.plot(t,models.cubic(*df_solution_cubic.values[0][:],t,y[-1],days[-1]))
        plt.plot(outliers[0,0,:],outliers[0,1,:],'ro',mfc='none',ms=10)
        plt.savefig("images/single_covid_cubic.pdf",bbox_inches = "tight")

    elif opt == 2:
        plt.plot(t,models.logistic(*df_solution_logistic.values[0][:],t))
        plt.plot(outliers[1,0,:],outliers[1,1,:],'ro',mfc='none',ms=10)
        plt.savefig("images/single_covid_logistic.pdf",bbox_inches = "tight")
    
    elif opt == 3:
        plt.plot(t,models.der_logistic(*df_solution_der_logistic.values[0][:],t))
        plt.plot(outliers[2,0,:],outliers[2,1,:],'ro',mfc='none',ms=10)
        plt.savefig("images/single_covid_der_logistic.pdf",bbox_inches = "tight")

    else:
        plt.plot(t,models.cubic(*df_solution_cubic.values[0][:],t,y[-1],days[-1]))
        plt.plot(t,models.logistic(*df_solution_logistic.values[0][:],t))
        plt.savefig("images/covid_all_models.pdf",bbox_inches = "tight")

    plt.show()
    plt.close()

# Reading the necessary data
df_solution_cubic       = pd.read_table("output/solutions_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_solution_logistic    = pd.read_table("output/solutions_covid_logistic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_solution_der_logistic= pd.read_table("output/solutions_covid_der_logistic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_train_set            = pd.read_table("data/covid_train.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_test_set             = pd.read_table("data/covid_test.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
n_train                 = int(df_train_set.values[0][0])
n_test                  = int(df_test_set.values[0][0])

with open("output/inf_sup_covid.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

# Numbers of previous days considered in model fitting
previous_days = int(xdata[0])

# Arrays allocation
y               = np.zeros(n_train)
y_later         = np.zeros(n_test)
days            = np.linspace(1,previous_days,previous_days)
days_later      = np.linspace(1 + previous_days,previous_days + n_test,n_test)
t               = np.linspace(1,previous_days + n_test,1000)

# Observation in the previous days considered
for i in range(n_train):
    y[i] = df_train_set.values[i+1][0]

# Observation in the following days
for i in range(n_test):
    y_later[i] = df_test_set.values[i+1][0]

with open("output/outliers_covid_cubic.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

noutliers       = int(xdata[0])
outliers        = np.empty((3,2,noutliers))
ind_outliers    = np.empty(noutliers,dtype=int)

for i in range(noutliers):
    ind_outliers[i] = int(xdata[i+1])

for i in range(noutliers):
    outliers[0,0,i] = days[ind_outliers[i]-1]
    outliers[0,1,i] = y[n_train - previous_days + ind_outliers[i]-1]

with open("output/outliers_covid_logistic.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

for i in range(noutliers):
    ind_outliers[i] = int(xdata[i+1])

for i in range(noutliers):
    outliers[1,0,i] = days[ind_outliers[i]-1]
    outliers[1,1,i] = y[n_train - previous_days + ind_outliers[i]-1]

# with open("output/outliers_covid_der_logistic.txt") as f:
#     lines = f.readlines()
#     xdata = [line.split()[0] for line in lines]

# for i in range(noutliers):
#     ind_outliers[i] = int(xdata[i+1])

# for i in range(noutliers):
#     outliers[2,0,i] = days[ind_outliers[i]-1]
#     outliers[2,1,i] = y[n_train - previous_days + ind_outliers[i]-1]

# plot_models(1)
plot_models(2)
# plot_models()
