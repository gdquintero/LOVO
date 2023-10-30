import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import models

def plot_models(opt=None):
    plt.rcParams.update({'font.size': 14})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlim(-3,max(later_days)+0.5)
    plt.ylim(min(min(y[len(y)-sup:]),min(y_later))-0.5,max(max(y[len(y)-sup:]),max(y_later))+0.5)
    
    if opt == 1:
        for i in range(sup - inf + 1):
            t = np.linspace(n-i,sup + n_test,h)
            plt.plot(t,models.cubic(*df_solution_cubic.values[i][:],t,y[-1],previous_days[-1]),label="days = " + str(inf+i))

    elif opt == 2:
        for i in range(sup - inf + 1):
            t = np.linspace(n-i,sup + n_test,h)
            plt.plot(t,models.logistic(*df_solution_logistic.values[i][:],t-n+i+1),label="days = " + str(inf+i))


    plt.plot(np.linspace(1,sup,sup),y[len(y)-sup:],"ko")
    plt.plot(later_days,y_later,"o",color="grey")
    plt.legend()
    plt.show()
    plt.close()

df_solution_cubic       = pd.read_table("output/solutions_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_solution_logistic    = pd.read_table("output/solutions_covid_logistic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_train_set            = pd.read_table("data/covid_train.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_test_set             = pd.read_table("data/covid_test.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
n_train                 = int(df_train_set.values[0][0])
n_test                  = int(df_test_set.values[0][0])

with open("output/inf_sup_covid.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

inf = int(xdata[0])
sup = int(xdata[1])
h = 1000
n = sup - inf + 1

# Arrays allocation
y               = np.zeros(n_train)
y_later         = np.zeros(n_test)
previous_days   = np.linspace(1,sup,sup)
later_days      = np.linspace(sup + 1,sup + n_test,n_test)


# Observation in the previous days considered
for i in range(n_train):
    y[i] = df_train_set.values[i+1][0]

# Observation in the following days
for i in range(n_test):
    y_later[i] = df_test_set.values[i+1][0]

plot_models(1)
# plot_models(2)

