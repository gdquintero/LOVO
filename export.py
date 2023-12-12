import numpy as np
import pandas as pd

df_accuracy_cubic = pd.read_table("output/accuracy_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_accuracy_cubic.to_excel("output/accuracy_cubic.xlsx")

df_accuracy_logistic = pd.read_table("output/accuracy_covid_logistic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_accuracy_logistic.to_excel("output/accuracy_logistic.xlsx")

df_cubic_latex = pd.read_table("output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_cubic_latex.to_excel("output/solution_cubic.xlsx")

df_cubic_latex = pd.read_table("data/cubic_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_cubic_latex.to_excel("data/cubic_latex.xlsx")

