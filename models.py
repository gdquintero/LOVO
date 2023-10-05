import numpy as np

# Farrington models
l = lambda t, a, b, c: (a * t - c) * np.exp(-b * t) + c

F = lambda t, a, b, c: 1.0 - np.exp((a/b) * t * np.exp(-b * t) + (1.0/b) * ((a/b) - c) * \
    (np.exp(-b * t) - 1.0) - c * t)

# Covid models
cubic = lambda x1,x2,x3,t,ym,tm: ym + x1 * (t - tm) + x2 * (t - tm)**2 + x3 * (t - tm)**3