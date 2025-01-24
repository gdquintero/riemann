import numpy as np
import scipy as sc
from tools import mobius_function

def H(alpha,t,k):
    s = alpha + t*1j
    return (1 - k**(1 - s)) * (sc.special.zeta(s) / s)

def Gn(alpha,t,n):
    res = 0
    for k in range(2,n):
        print(k)
        res += (mobius_function(k)/k) * H(alpha,t,k)

    return res + 1/(alpha + t*1j)

def func(x,y):
    return (x**2) * y

n = 100000
alpha = 2

res, _ =  sc.integrate.quad(lambda t: np.absolute(Gn(alpha,t,n))**2,-np.inf,np.inf)

print(res)
