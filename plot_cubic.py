import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_table("data/cubic.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)

plt.plot(df[0].values,df[1].values,"o")
plt.show()