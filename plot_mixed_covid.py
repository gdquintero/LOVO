import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import models

def plot_models(opt=None):
    plt.plot(previous_days,y[:sup],"ko")
    plt.plot(later_days,y_later,"o",color="grey")
    # plt.ylim(min(min(y[n_train - inf:]),min(y_later))-0.5,max(max(y[n_train - inf:]),max(y_later))+0.5)

    if opt == 1:
        for i in range(sup-inf):
            plt.plot(t,models.cubic(*df_solution_cubic.values[i][:],t,y[-1],previous_days[-1]))

        plt.savefig("images/single_covid_cubic.pdf",bbox_inches = "tight")

    elif opt == 2:
        for i in range(sup-inf):
            plt.plot(t,models.logistic(*df_solution_logistic.values[i][:],t))

        plt.savefig("images/single_covid_logistic.pdf",bbox_inches = "tight")

    else:
        plt.plot(t,models.cubic(*x[0,:],t,y[-1],previous_days[-1]))
        plt.plot(t,models.logistic(*x[1,:],t))
        plt.savefig("images/covid_all_models.pdf",bbox_inches = "tight")

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

# Arrays allocation
y               = np.zeros(n_train)
y_later         = np.zeros(n_test)
previous_days   = np.linspace(1,sup,sup)
later_days      = np.linspace(sup + 1,sup + n_test,n_test)
t               = np.linspace(1,sup + n_test,1000)

# Observation in the previous days considered
for i in range(n_train):
    y[i] = df_train_set.values[i+1][0]

# Observation in the following days
for i in range(n_test):
    y_later[i] = df_test_set.values[i+1][0]



plot_models(2)

