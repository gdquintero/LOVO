import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

plt.rcParams.update({'font.size': 14})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
NUM_COLORS = 11
cm = plt.get_cmap('rainbow')
cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle(color=[scalarMap.to_rgba(NUM_COLORS-i) for i in range(NUM_COLORS)])

def poly(x,t):
    return x[0] + x[1] * t + x[2] * (t**2) + x[3] * (t**3)

df_data = pd.read_table("~/github/LOVO/data/cubic.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)
df_sol = pd.read_table("~/github/LOVO/output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)

plt.plot(df_data[0].values,df_data[1].values,"ok")

# for i in range(len(df_sol)):
#     plt.plot(t,poly(df_sol.values[i][:4],t),label="o = " + str(i))


# plt.plot(t,poly(df_sol.values[0][:4],t),label="o = 0")
# plt.plot(t,poly(df_sol.values[5][:4],t),label="o = 5")
plt.legend()
plt.show()
