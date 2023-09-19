import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df_solution = pd.read_table("output/solution_covid_cubic.txt",delimiter=" ",header=None,skiprows=0)

print(df_solution)