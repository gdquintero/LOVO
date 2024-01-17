import numpy as np
import pandas as pd

# df = pd.read_table("output/accuracy_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
# df.to_excel("output/accuracy_cubic.xlsx")

# df = pd.read_table("output/accuracy_covid_logistic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
# df.to_excel("output/accuracy_logistic.xlsx")

# df= pd.read_table("output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
# df.to_excel("output/solution_cubic.xlsx")

# df = pd.read_table("data/cubic_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
# df.to_excel("data/cubic_latex.xlsx")

df= pd.read_table("output/measles_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel("output/measles_latex.xlsx")

df= pd.read_table("output/mumps_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel("output/mumps_latex.xlsx")

df= pd.read_table("output/rubella_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel("output/rubella_latex.xlsx")