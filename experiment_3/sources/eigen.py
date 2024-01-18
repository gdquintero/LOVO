import numpy as np

# Definir el valor de t
t = 3  # Puedes cambiar este valor

# Definir la matriz A
A = np.array([[1, t, t**2, t**3],
              [t, t**2, t**3, t**4],
              [t**2, t**3, t**4, t**5],
              [t**3, t**4, t**5, t**6]])

# Calcular autovalores
autovalores = np.linalg.eigvals(A)

print("Autovalores:", autovalores)