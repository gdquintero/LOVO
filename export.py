import numpy as np
import pandas as pd

df_accuracy_cubic = pd.read_table("output/accuracy_covid_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

df_accuracy_cubic.to_excel("output/accuracy.xlsx")