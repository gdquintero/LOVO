import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def poly(x,t):
    return x[0] + x[1] * t + x[2] * (t**2) + x[3] * (t**3)

df_data = pd.read_table(parent+"/data/cubic.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_cubic.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

with open(parent+"/output/num_mixed_test.txt") as f:
    lines = f.readlines()
    lim = [line.split()[0] for line in lines]

inf = int(lim[0])
sup = int(lim[1])
n = sup - inf + 1

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
NUM_COLORS = n
cm = plt.get_cmap('rainbow')
cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS)
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle(color=[scalarMap.to_rgba(NUM_COLORS-i) for i in range(NUM_COLORS)])

t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)

for i in range(n):
    plt.plot(t,poly(df_sol.values[i][:4],t),label="o = "+str(inf+i),linewidth=1.5)
    plt.legend(loc='upper left')

l = plt.plot(df_data[0].values,df_data[1].values,"ok")
plt.setp(l, 'markersize', 4)
plt.xticks(range(-1, 4, 1))
plt.yticks(range(-4, 5, 2))
plt.show()
