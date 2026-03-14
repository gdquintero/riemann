import scipy as sc
import time

def H(alpha,t,k):
    s = alpha + t*1j
    return (1 - k**(1 - s)) * (sc.special.zeta(s) / s)


print(sc.special.zeta(0.5+100000000j))


