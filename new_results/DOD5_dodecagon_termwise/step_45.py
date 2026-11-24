import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
mp.mp.dps = 700   # headroom
xc = pickle.load(open('dod5_xi36_closed.pkl','rb'))
Rc = pickle.load(open('dod5_R_content.pkl','rb'))
xv = pickle.load(open(D2+"Xi_verified_620.pkl","rb"))

# --- Xi36 part value (recompute at 700 for headroom) ---
C = consts(700)
LI3B = [Li3v(z) for z in BASIS3_ZS]
PBV = [prodval(m, C) for m in PROD_META]
IBV = [C['z3'], C['L3c'], C['z2']*C['l2'], C['z2']*C['l3'], C['z2']*C['lu']]
def fval(p):
    return mp.mpf(p[0].numerator)/p[0].denominator + mp.sqrt(3)*mp.mpf(p[1].numerator)/p[1].denominator
xi36v = mp.mpf(0)
for j in range(8): xi36v += fval(xc['A3'][j])*LI3B[j]
for j in range(19): xi36v += fval(xc['AP'][j])*PBV[j]
for m in range(5): xi36v += fval(xc['AI'][m])*IBV[m]

# --- R part value ---
import sympy as sp
l2s,l3s,lus = sp.symbols('l2 l3 lu')
Li2A = sp.Function('Li2A'); Li3A = sp.Function('Li3A')
r3m = mp.sqrt(3)
vals_sym = {l2s: mp.log(2), l3s: mp.log(3), lus: mp.log(2-r3m)}  # dod3_R convention
def Li3Rv(z): return mp.re(mp.polylog(3, mp.mpc(mpf_of(z))))
Rv3 = mp.mpf(0)
for pt,c in Rc['li3_content'].items():
    Rv3 += fval(psplit(c))*Li3Rv(pt)
for s in ('li2log','logcube'):
    e = sp.sympify(Rc[s])
    e = e.subs(vals_sym)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,700))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,700)))))
    Rv3 += mp.mpf(sp.N(e, 700))
print("Xi36 closed value:", mp.nstr(xi36v,30))
print("R sum value      :", mp.nstr(Rv3,30))
Xi_closed = xi36v + Rv3
print("Xi closed        :", mp.nstr(Xi_closed,40))
print("Xi_from_S        :", mp.nstr(mp.mpf(xv['Xi_from_S']),40))
print("|diff| =", mp.nstr(abs(Xi_closed - mp.mpf(xv['Xi_from_S'])),3), "  t=%.0fs"%(time.time()-t0))

# --- P'''(1) ---
L = mp.sqrt(26-15*r3m)
P3 = mp.mpf(pickle.load(open(D1+"P3_630.pkl","rb")))
S_nod2 = mp.mpf(xv['S_nod2'])
P3_closed = 72*L*(S_nod2 + Xi_closed)
print("\nL =", mp.nstr(L,25))
print("P3 closed =", mp.nstr(P3_closed, 40))
print("P3_630    =", mp.nstr(P3, 40))
print("|P3_closed - P3_630| =", mp.nstr(abs(P3_closed-P3),3))
pickle.dump({'xi36v':xi36v,'Rv3':Rv3,'Xi_closed':Xi_closed,'P3_closed':P3_closed}, open('dod5_xi_p3_closed.pkl','wb'))
print("saved, elapsed %.0fs"%(time.time()-t0))
