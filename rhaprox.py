import numpy as np
import scipy as sc
from tools import mobius_function

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

values = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]


for n in values:
    res = 0
    res, _ =  sc.integrate.quad(lambda t: np.absolute(Gn(alpha,t,n))**2,-np.inf,np.inf)

    print("n = " + str(n) + " | ","Integral =", str(res))
