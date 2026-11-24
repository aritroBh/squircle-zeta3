import sys, time, itertools, pickle
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values
TW = TW3()
Gsym = {g: col_Li3(sp.Rational(g), SS235) for g in Gargs}
PBspec = [('-1/3',2),('-1/3',3),('-1/3',5),('-1/4',2),('-1/4',3),('-1/4',5),('1/4',2),('1/4',3),('1/4',5),
          ('-5/3',2),('-5/3',3),('-5/3',5),('-5/4',2),('-5/4',3),('-5/4',5),('1/5',3),('1/5',5),('5/8',2),('5/8',3)]
PB = [(sp.Rational(w), p) for w,p in PBspec]
Pmat = sp.Matrix.hstack(*[col_prod(w, p, SS235) for w,p in PB])
sel = ['1/3', '2/3', '2/5', '-1/5', '1/10', '-5/3', '5', '3/5']
selq = [sp.Rational(g) for g in sel]
M = sp.Matrix.hstack(*[Gsym[g] for g in selq], Pmat)
sols = {e: [sp.Rational(c) for c in M.LUsolve(sp.Matrix(TW[e]))] for e in (1,2,3)}

DPS = 300
mp.mp.dps = DPS
def Li3v(a): return mp.re(mp.polylog(3, mp.mpf(a)))
def Li2v(a): return mp.re(mp.polylog(2, mp.mpf(a)))
z2 = mp.pi**2/6; z3 = mp.zeta(3)
lam = pickle.load(open(base+"OCT9c_octagon_close/lam_mine.pkl","rb"))
for e in (1,2,3):
    eng = get_engine(e)
    NUM = atom_values(eng, e, lam, DPS)
    W3v = mp.mpf(sp.N(D['s123'][e]['W3'].subs(list(NUM.items())), DPS))
    c = sols[e]
    li = mp.fsum([mp.mpf(str(c[i]))*Li3v(selq[i]) for i in range(8)])
    pr = mp.fsum([mp.mpf(str(c[8+j]))*(-Li2v(PB[j][0])*mp.log(PB[j][1])) for j in range(19)])
    R = W3v - li - pr
    vec = [R, z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)]
    rel = mp.pslq(vec, tol=mp.mpf('1e-250'), maxcoeff=10**30, maxsteps=10000)
    res = abs(mp.fsum([mp.mpf(r)*v for r,v in zip(rel,vec)])) if rel else None
    print(f"e{e}: R={mp.nstr(R,20)}")
    print(f"   rel={rel}")
    print(f"   residual={mp.nstr(res,2) if res else 'NONE'}")
