import sys, os, time
work = '<project>/research/08_theory/TOW2R_review/tow2_copy'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
i = T.i

def feed(delta):
    tot = 0
    for dp in range(0, delta+1):
        cA = T.tau1(delta-dp) * (m+1) / (2*d-dp-1)
        cB = T.tau2(delta-dp) * (d-m+1) / (2*d-dp-1)
        tot += cA * T.arho(dp).subs({d: d+1, m: m+1}) * T.rA0()
        tot += cB * T.arho(dp).subs({d: d+1, m: m}) * T.rB0()
    return tot

t0 = time.time()
rho4 = T.arho(4)
g = d**3 - 9*d**2 + sp.Rational(113,4)*d - i*d/2 + i*m - 29
e41 = g/((2*d-3)*(2*d-5))
lhs4 = rho4 + e41
print('built lhs in', round(time.time()-t0,1),'s; starting cancel...')
rem4 = sp.cancel(sp.together(lhs4 - feed(4)))
print('V(4) non-adjacent remainder:', rem4, ' (', round(time.time()-t0,1),'s )')
