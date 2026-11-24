import sys, time, itertools, pickle, random
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values
t0=time.time()
TW = TW3()
Gsym = {g: col_Li3(sp.Rational(g), SS235) for g in Gargs}
PB = [(sp.Rational(w),p) for w,p in [('-1/3',2),('-1/3',3),('-1/3',5),('-1/4',2),('-1/4',3),('-1/4',5),('1/4',2),('1/4',3),('1/4',5),
      ('-5/3',2),('-5/3',3),('-5/3',5),('-5/4',2),('-5/4',3),('-5/4',5),('1/5',3),('1/5',5),('5/8',2),('5/8',3)]]
PMAT = sp.Matrix.hstack(*[col_prod(w,p,SS235) for w,p in PB])
DPS = 250
mp.mp.dps = DPS
def Li3v(a): return mp.re(mp.polylog(3, mp.mpf(a)))
def Li2v(a): return mp.re(mp.polylog(2, mp.mpf(a)))
lam = pickle.load(open(base+"OCT9c_octagon_close/lam_mine.pkl","rb"))
W3v = {}
for e in (1,2,3):
    eng = get_engine(e)
    NUM = atom_values(eng, e, lam, DPS)
    W3v[e] = mp.mpf(sp.N(D['s123'][e]['W3'].subs(list(NUM.items())), DPS))
pvals = [(-Li2v(w))*mp.log(p) for w,p in PB]
z2 = mp.pi**2/6; z3 = mp.zeta(3)
def heights_G8(sel):
    selq = [sp.Rational(g) for g in sel]
    M = sp.Matrix.hstack(*[Gsym[g] for g in selq], PMAT)
    if M.rank() < 27: return None
    lvals = [Li3v(g) for g in selq]
    out = []
    for e in (1,2,3):
        x = M.LUsolve(sp.Matrix(TW[e]))
        if M*x != sp.Matrix(TW[e]): return None
        c = [mp.mpf(str(sp.Rational(cc))) for cc in x]
        R = W3v[e] - mp.fsum([c[i]*lvals[i] for i in range(8)]) - mp.fsum([c[8+j]*pvals[j] for j in range(19)])
        vec = [R, z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)]
        rel = mp.pslq(vec, tol=mp.mpf('1e-200'), maxcoeff=10**30, maxsteps=10000)
        if not rel: return None
        out.append((R, rel, abs(mp.fsum([mp.mpf(r)*v for r,v in zip(rel,vec)]))))
    return out
bases = [
 ['1/3','2/3','2/5','-1/5','1/10','-5/3','5','3/5'],
 ['10/9','-2/3','-1/2','1/3','-5/3','-1/5','5/4','2/5'],
 ['1/3','2/5','5/8','-5/3','2/3','9/5','-2/3','4/5'],
 ['-2/3','5/6','-1/2','-5/3','-5/4','1/5','5/8','-1/9'],
]
for sel in bases:
    h = heights_G8(sel)
    if h is None: print(sel, "FAIL"); continue
    mx = max(max(abs(r) for r in x[1]) for x in h)
    print(sel, "->", [max(abs(r) for r in x[1]) for x in h], "max", f"{mx:.2e}", f"({time.time()-t0:.1f}s)")
