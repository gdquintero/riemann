import scipy as sc

def H(alpha,t,k):
    s = alpha + t*1j
    return (1 - k**(1 - s)) * (sc.special.zeta(s) / s)


alpha = 2
t = 3
k = 100

print(H(alpha,t,k))

