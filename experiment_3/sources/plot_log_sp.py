import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def plot1():
    fig,ax = plt.subplots(1, 1)
    ax.tick_params(axis='both',direction='in',which='both')
    ax.loglog(x,y,":o",color="darkgreen",lw=1,ms=4)
    ax.set_xscale('linear')
    plt.ylim(0,10)
    # ax.set_ylim(min(np.log10(y)),max(np.log10(y))+10)
    plt.xlabel("Number of outliers $o$")
    plt.ylabel("$S_q(x^*)$ (log scale)")
    plt.savefig(parent+"/images/log10.pdf",bbox_inches="tight")
    
    # plt.show()

def plot2():
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    ax.semilogx(x,y)
    ax.set_yscale('linear')
    plt.show()

df_sp = pd.read_table(parent+"/output/log_sp.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

with open(parent+"/output/num_mixed_test.txt") as f:
    lines = f.readlines()
    lim = [line.split()[0] for line in lines]

inf = int(lim[0])
sup = int(lim[1])
n = sup - inf + 1
x = np.linspace(inf,sup,n)
y = df_sp[0].values

plot1()
