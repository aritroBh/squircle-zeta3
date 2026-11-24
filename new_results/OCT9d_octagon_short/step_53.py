import sys, time, itertools, pickle
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values
t0=time.time()
basis0 = [sp.Rational(g) for g in pickle.load(open("basis0.pkl","rb"))['basis0']]
PB0 = [(sp.Rational(w),p) for w,p in pickle.load(open("PB0.pkl","rb"))['PB0']]
sol0 = [sp.Rational(c) for c in pickle.load(open("sol0_short.pkl","rb"))['sol0']]
DPS = 250
mp.mp.dps = DPS
def Li3v(a): return mp.re(mp.polylog(3, mp.mpf(a)))
def Li2v(a): return mp.re(mp.polylog(2, mp.mpf(a)))
z2 = mp.pi**2/6; z3 = mp.zeta(3)
lam = pickle.load(open(base+"OCT9c_octagon_close/lam_mine.pkl","rb"))
eng = get_engine(0)
NUM = atom_values(eng, 0, lam, DPS)
W3v = mp.mpf(sp.N(D['Wp0']['W3_e0'].subs(list(NUM.items())), DPS))
li = mp.fsum([mp.mpf(str(sol0[i]))*Li3v(basis0[i]) for i in range(20)])
pr = mp.fsum([mp.mpf(str(sol0[20+j]))*(-Li2v(PB0[j][0])*mp.log(PB0[j][1])) for j in range(44)])
R = W3v - li - pr
print("R_0 =", mp.nstr(R, 30))
vec = [R, z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5), z2*mp.log(7)]
rel = mp.pslq(vec, tol=mp.mpf('1e-200'), maxcoeff=10**30, maxsteps=10000)
print("rel:", rel)
if rel:
    print("residual:", mp.nstr(abs(mp.fsum([mp.mpf(r)*v for r,v in zip(rel,vec)])), 2))
print(f"({time.time()-t0:.1f}s)")
