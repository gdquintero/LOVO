import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import models
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

def plot_mixed(ind,t,inf,df_seropositives,df_mixed):
    size_img = 0.6
    plt.rcParams.update({'font.size': 11})
    plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    NUM_COLORS = 11
    cm = plt.get_cmap('rainbow')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS+1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(color=[scalarMap.to_rgba(NUM_COLORS-i) for i in range(NUM_COLORS)])
    plt.ylim([0,1.05])
    

    for i in range(n):
        plt.plot(t,models.F(t,*df_mixed.iloc[i].values),label="o = "+str(inf+i),linewidth=1)
        plt.legend(loc='lower right',prop={'size':7})

    plt.plot(df_seropositives[0].values,df_seropositives[ind].values,"ko",ms=3)
    plt.tick_params(axis='both',direction='in')
    plt.xticks(range(0, 71, 10))
    # plt.yticks(range(-4, 5, 2))
    plt.savefig(parent+"/images/"+str(ind)+".pdf",bbox_inches = "tight")
    
    # plt.show()
    plt.close()

df_seropositives = pd.read_table(parent+"/data/seropositives.txt",delimiter=" ",header=None,skiprows=1)
df_mixed_measles = pd.read_table(parent+"/output/solutions_mixed_measles.txt",delimiter=" ",header=None,skiprows=0)
df_mixed_mumps   = pd.read_table(parent+"/output/solutions_mixed_mumps.txt",delimiter=" ",header=None,skiprows=0)
df_mixed_rubella = pd.read_table(parent+"/output/solutions_mixed_rubella.txt",delimiter=" ",header=None,skiprows=0)

with open(parent+"/output/num_mixed_test.txt") as f:
    lines = f.readlines()
    lim = [line.split()[0] for line in lines]

inf = int(lim[0])
sup = int(lim[1])
n = sup - inf + 1
t = np.linspace(0,70,1000)

plot_mixed(1,t,inf,df_seropositives,df_mixed_measles)
plot_mixed(2,t,inf,df_seropositives,df_mixed_mumps)
plot_mixed(3,t,inf,df_seropositives,df_mixed_rubella)


