import numpy as np
import scipy as sc
from prime_tools import mobius_function

def H(alpha,t,k):
    s = alpha + t*1j
    return (1 - k**(1 - s)) * (sc.special.zeta(s) / s)

def Gn(alpha,t,n):
    res = 0
    for k in range(2,n):
        res += (mobius_function(k)/k) * H(alpha,t,k)

    return res + 1/(alpha + t*1j)

def func(x,y):
    return (x**2) * y

alpha = 0.9

values = [1000]


for n in values:
    res, err = 0, 0
    res, err =  sc.integrate.quad(lambda t: np.absolute(Gn(alpha,t,n))**2,-np.inf,np.inf,epsrel=1e-12, epsabs=1e-15)

print(res,err)