from mpmath import mp, mpf, pi, sin, cos, quad
import time, json, os
mp.dps = 45

def pi_p(p):
    p = mpf(p)
    return 2*pi/(p*sin(pi/p))

def F_p(y, p):
    p=mpf(p); y=mpf(y)
    if y == 0: return mpf(0)
    if y == 1: return pi_p(p)/2
    return y*quad(lambda s: (1-(y*s)**p)**(-1/p), [0,1])

def Phi_p(y,p):
    return pi*F_p(y,p)/pi_p(p)

def b_n_quad(n, p):
    return 4/(n*pi)*quad(lambda y: cos(n*Phi_p(y,p)), [0,1])

# quick timing check at p near 2
t0=time.time(); print(mp.nstr(b_n_quad(3, 2.02),30), f"{time.time()-t0:.1f}s")
