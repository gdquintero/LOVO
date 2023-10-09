import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import models

df_solutions = pd.read_table("output/solutions_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_train_set = pd.read_table("data/covid_train.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_test_set = pd.read_table("data/covid_test.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
n_train = int(df_train_set.values[0][0])
n_test = int(df_test_set.values[0][0])

with open("output/inf_sup_covid.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

inf = int(xdata[0])
sup = int(xdata[1])

n = sup - inf + 1

x   = np.zeros(3)
y   = np.zeros(n_train)
y_later = np.zeros(n_test)
days = np.linspace(1,inf,inf)
days_later = np.linspace(1 + inf,inf + n_test,n_test)
total_days = np.linspace(1,sup - inf + n_test+1,sup - inf + n_test+1)
t   = np.linspace(1,sup - inf + n_test,1000)


for i in range(n_train):
    y[i] = df_train_set.values[i+1][0]

for i in range(n_test):
    y_later[i] = df_test_set.values[i+1][0]


print(n)