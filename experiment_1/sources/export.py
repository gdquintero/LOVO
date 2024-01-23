import numpy as np
import pandas as pd
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

df= pd.read_table(parent+"/output/measles_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/measles_latex.xlsx")

df= pd.read_table(parent+"/output/mumps_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/mumps_latex.xlsx")

df= pd.read_table(parent+"/output/rubella_latex.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/rubella_latex.xlsx")

df= pd.read_table(parent+"/output/measles_latex2.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/measles_latex2.xlsx")

df= pd.read_table(parent+"/output/mumps_latex2.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/mumps_latex2.xlsx")

df= pd.read_table(parent+"/output/rubella_latex2.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df.to_excel(parent+"/output/rubella_latex2.xlsx")