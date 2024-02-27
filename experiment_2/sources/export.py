import numpy as np
import pandas as pd
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df= pd.read_table(parent+"/output/latex_5.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/latex_5.xlsx")

df= pd.read_table(parent+"/output/latex_10.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/latex_10.xlsx")

df= pd.read_table(parent+"/output/latex_15.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/latex_15.xlsx")

df= pd.read_table(parent+"/output/latex_20.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/latex_20.xlsx")

df= pd.read_table(parent+"/output/latex_25.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/latex_25.xlsx")

df= pd.read_table(parent+"/output/latex_30.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/latex_30.xlsx")

df= pd.read_table(parent+"/output/mean_per_day_re.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/mean_per_day_re.xlsx")

df= pd.read_table(parent+"/output/successfull_threshold.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/successfull_threshold.xlsx")