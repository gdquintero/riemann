import numpy as np
import scipy as sc
import time
from prime_tools import mobius_function

def H(alpha,t,k):
    s = alpha + t*1j
    return (1 - k**(1 - s)) * (sc.special.zeta(s) / s)

def Gn(alpha, t, n):
    # Precalcular valores de mobius_function(k) / k
    k_vals = np.arange(2, n)
    mobius_vals = np.array([mobius_function(k) for k in k_vals]) / k_vals

    # Vectorizar cálculo de H(alpha, t, k)
    H_vals = np.array([H(alpha, t, k) for k in k_vals])

    # Calcular la suma
    res = np.sum(mobius_vals * H_vals)

    # Añadir el término final
    return res + 1 / (alpha + t * 1j)

# Acumular resultados para escribir en un solo paso
results = []

param = float(np.loadtxt("parameters.txt"))

# tol = 5e-2
# max_iter = 1000

# # Tiempo máximo 10 min
# max_time = 60 * 10

# points = 100

# # Generamos particion del intervalo cerrado [0.5,1]
# partition = np.linspace(0.5,1,points)

# for alpha in partition:

#     start = time.time() # Marca el tiempo de inicio
#     iter = 0
#     n = 0
#     res = np.inf

#     while True:
#         if time.time() - start > max_time:
#             break
#         elif iter > max_iter:
#             break
#         elif res <= tol:
#             break
#         n += 1
#         iter += 1
        
#         # Realizar la integración con scipy
#         res, err = sc.integrate.quad(
#             lambda t: np.absolute(Gn(alpha, t, n)) ** 2,
#             -np.inf,
#             np.inf,
#             epsrel=1e-8,
#             epsabs=1e-8,
#         )  

#     finish = time.time() # Marca el tiempo de finalizacion
#     total_time = finish - start
#     # Acumular resultado
#     results.append(f"{alpha} {res} {n} {total_time}\n")

#     # Escribir todo al archivo en un solo paso
#     with open("output_tol_0.05.txt", "w") as f:
#         f.writelines(results)

#     start, finish = 0, 0

