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

alpha = float(np.loadtxt("alpha.txt"))

for n in range(36000, 41000, 1000):

    start = time.time() # Marca el tiempo de inicio
    # Realizar la integración con scipy
    res, err = sc.integrate.quad(
        lambda t: np.absolute(Gn(alpha, t, n)) ** 2,
        -np.inf,
        np.inf,
        epsrel=1e-8,
        epsabs=1e-8,
    )

    finish = time.time() # Marca el tiempo de finalizacion
    total_time = finish - start
    # Acumular resultado
    results.append(f"{n} {alpha} {res} {total_time}\n")

# Escribir todo al archivo en un solo paso
with open("output8.txt", "w") as f:
    f.writelines(results)