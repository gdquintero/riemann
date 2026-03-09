import scipy as sc
import time

def H(alpha,t,k):
    s = alpha + t*1j
    return (1 - k**(1 - s)) * (sc.special.zeta(s) / s)


alpha = 0.9
t = 2738947983427.0
n = 2


print(sc.special.zeta(1+1j))


