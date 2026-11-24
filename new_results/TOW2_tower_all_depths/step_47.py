import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
i = T.i

t0 = time.time()
e41 = T.e_desc(4, 1)
print('e_desc(4,1) computed in', round(time.time()-t0,1),'s')
# compare with fitted form
g = d**3 - 9*d**2 + sp.Rational(113,4)*d - i*d/2 + i*m - 29
den = (2*d-5)*(2*d-3)
fit = g/den
print('e_desc(4,1) == fitted e(4,1):', sp.cancel(sp.together(e41 - fit)) == 0)

e30 = T.e_desc(3, 0)
rho3 = T.arho(3)
print('V(3): e_desc(3,0) == rho(3):', sp.cancel(sp.together(e30 - rho3)) == 0)

e40 = T.e_desc(4, 0)
rho4 = T.arho(4)
print('V(4): e_desc(4,0) == rho(4):', sp.cancel(sp.together(e40 - rho4)) == 0)
print('total time', round(time.time()-t0,1))
