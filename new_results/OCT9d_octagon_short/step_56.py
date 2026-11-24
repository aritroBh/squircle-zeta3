import sys, time, pickle
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values
t0=time.time()
SF = pickle.load(open("short_form.pkl","rb"))
G8 = [sp.Rational(g) for g in SF['G8']]
PB19 = [(sp.Rational(w),p) for w,p in SF['PB19']]
B0 = [sp.Rational(g) for g in SF['B0']]
PB0 = [(sp.Rational(w),p) for w,p in SF['PB0']]
sols123 = {e: [sp.Rational(c) for c in SF['sols123'][e]] for e in (1,2,3)}
sol0 = [sp.Rational(c) for c in SF['sol0']]
lam_hi = pickle.load(open("lam_hi.pkl","rb"))

def Li3v(a): return mp.re(mp.polylog(3, mp.mpf(a)))
def Li2v(a): return mp.re(mp.polylog(2, mp.mpf(a)))

def compute_R(e, dps):
    mp.mp.dps = dps
    eng = get_engine(e)
    NUM = atom_values(eng, e, lam_hi, dps)
    if e == 0:
        W3v = mp.mpf(sp.N(D['Wp0']['W3_e0'].subs(list(NUM.items())), dps))
        c = sol0
        li = mp.fsum([mp.mpf(str(c[i]))*Li3v(B0[i]) for i in range(20)])
        pr = mp.fsum([mp.mpf(str(c[20+j]))*(-Li2v(PB0[j][0])*mp.log(PB0[j][1])) for j in range(44)])
    else:
        W3v = mp.mpf(sp.N(D['s123'][e]['W3'].subs(list(NUM.items())), dps))
        c = sols123[e]
        li = mp.fsum([mp.mpf(str(c[i]))*Li3v(G8[i]) for i in range(8)])
        pr = mp.fsum([mp.mpf(str(c[8+j]))*(-Li2v(PB19[j][0])*mp.log(PB19[j][1])) for j in range(19)])
    return W3v - li - pr

def pslq_const(R, e, dps):
    mp.mp.dps = dps
    z2 = mp.pi**2/6; z3 = mp.zeta(3)
    vec = [R, z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)] + ([z2*mp.log(7)] if e==0 else [])
    rel = mp.pslq(vec, tol=mp.mpf(10)**(-(dps-90)), maxcoeff=10**30, maxsteps=20000)
    res = abs(mp.fsum([mp.mpf(r)*v for r,v in zip(rel,vec)])) if rel else None
    return rel, res

results = {}
for e in (0,1,2,3):
    Rlo = compute_R(e, 310); rello, reslo = pslq_const(Rlo, e, 310)
    Rhi = compute_R(e, 410); relhi, reshi = pslq_const(Rhi, e, 410)
    same = (rello == relhi)
    results[e] = (Rhi, relhi, reslo, reshi, same)
    print(f"edge {e}: identical at 310 & 410: {same}")
    print(f"  R = {mp.nstr(Rhi, 25)}")
    print(f"  rel = {relhi}")
    print(f"  residual 310: {mp.nstr(reslo,2)}, residual 410: {mp.nstr(reshi,2)}  ({time.time()-t0:.1f}s)")
pickle.dump({e: (mp.nstr(results[e][0], 420), results[e][1]) for e in results}, open("invis_short.pkl","wb"))
print("saved invis_short.pkl")
