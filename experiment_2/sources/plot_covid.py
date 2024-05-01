import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as mplcm
import matplotlib.colors as colors

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

plt.rcParams.update({'font.size': 13})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')