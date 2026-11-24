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
EP = set()
for e in (1,2,3):
    for z in edge_points(e): EP.add(z)
ORB = set()
for z in list(EP): ORB |= orbit(z)
ORBSM = set()
for z in list(SMALL235): ORBSM |= orbit(z)
def nice(g):
    g = Fr(g)
    if g in EP: return 0
    if g in ORB: return 1
    if g in SMALL235: return 2
    if g in ORBSM: return 3
    return 4
wlist = sorted([g for g in Gargs], key=lambda g:(nice(g), max(abs(Fr(g).numerator),abs(Fr(g).denominator)), str(g)))
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
l3cache = {}
def Li3v_c(g):
    if g not in l3cache: l3cache[g] = Li3v(sp.Rational(g))
    return l3cache[g]
def heights_G8(sel):
    selq = [sp.Rational(g) for g in sel]
    M = sp.Matrix.hstack(*[Gsym[g] for g in selq], PMAT)
    if M.rank() < 27: return None
    lvals = [Li3v_c(g) for g in selq]
    hs = []
    for e in (1,2,3):
        x = M.LUsolve(sp.Matrix(TW[e]))
        if M*x != sp.Matrix(TW[e]): return None
        c = [mp.mpf(str(sp.Rational(cc))) for cc in x]
        R = W3v[e] - mp.fsum([c[i]*lvals[i] for i in range(8)]) - mp.fsum([c[8+j]*pvals[j] for j in range(19)])
        vec = [R, z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)]
        rel = mp.pslq(vec, tol=mp.mpf('1e-200'), maxcoeff=10**30, maxsteps=10000)
        if not rel: return None
        hs.append(max(abs(r) for r in rel))
    return hs
sel = ['10/9','-2/3','-1/2','1/3','-5/3','-1/5','5/4','2/5']
best = heights_G8(sel)
bestmx = max(best)
print("start:", best, f"{bestmx:.2e}")
improved = True
it = 0
while improved and it < 2:
    improved = False; it += 1
    for i in range(8):
        for g in wlist:
            g = str(g)
            if g in sel: continue
            sel2 = sel[:i]+[g]+sel[i+1:]
            h = heights_G8(sel2)
            if h is None: continue
            if max(h) < bestmx:
                sel, best, bestmx = sel2, h, max(h); improved = True
        print(f"pass {it} i={i}: {['%.1e'%x for x in best]}  sel={sel} ({time.time()-t0:.1f}s)")
print("FINAL sel:", sel, "heights:", best)
pickle.dump({'G8_final': sel}, open("G8_final.pkl","wb"))
