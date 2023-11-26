import numpy as np
import pandas as pd

df = pd.read_table("output/mean_percentage_error.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

print(sum(df.values / 100.))