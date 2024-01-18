import numpy as np

# Farrington models
def l(t,a,b,c): 
    return (a * t - c) * np.exp(-b * t) + c

def F(t,a,b,c):
    return 1.0 - np.exp((a/b) * t * np.exp(-b * t) + \
    (1.0/b) * ((a/b) - c) * (np.exp(-b * t) - 1.0) - c * t)

# Covid models
def cubic(x1,x2,x3,t,ym,tm):
    return ym + x1 * (t - tm) + x2 * (t - tm)**2 + x3 * (t - tm)**3

    