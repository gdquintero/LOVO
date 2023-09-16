import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import models

def plot_mixed(ind,t,inf,df_seropositives,df_mixed):
    disease = [r"Measles",r"Mumps",r"Rubella"]
    plt.rcParams.update({'font.size': 12})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # plt.figure(figsize=(10,7))
    # plt.xticks(fontsize=18)
    # plt.yticks(fontsize=18)
    plt.ylim([0,1.05])
    

    for i in range(n):
        plt.plot(t,models.F(t,*df_mixed.iloc[i].values),label="o = "+str(inf+i),linewidth=1.5)
        # plt.title(disease[ind-1],fontsize = 18)
        # plt.title(disease[ind-1])
        plt.legend( loc='lower right')

    l = plt.plot(df_seropositives[0].values,df_seropositives[ind].values,"ko")
    plt.setp(l, 'markersize', 6)

    plt.savefig(disease[ind-1]+".pdf",bbox_inches = "tight")
    plt.show()
    plt.close()

df_seropositives = pd.read_table("output/seropositives.txt",delimiter=" ",header=None,skiprows=1)
df_mixed_measles = pd.read_table("output/solutions_mixed_measles.txt",delimiter=" ",header=None,skiprows=0)
df_mixed_mumps   = pd.read_table("output/solutions_mixed_mumps.txt",delimiter=" ",header=None,skiprows=0)
df_mixed_rubella = pd.read_table("output/solutions_mixed_rubella.txt",delimiter=" ",header=None,skiprows=0)

with open("output/num_mixed_test.txt") as f:
    lines = f.readlines()
    lim = [line.split()[0] for line in lines]

inf = int(lim[0])
sup = int(lim[1])
n = sup - inf + 1
t = np.linspace(0,70,1000)

plot_mixed(1,t,inf,df_seropositives,df_mixed_measles)
plot_mixed(2,t,inf,df_seropositives,df_mixed_mumps)
plot_mixed(3,t,inf,df_seropositives,df_mixed_rubella)


