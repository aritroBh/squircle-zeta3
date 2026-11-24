import sys, time, itertools, pickle, random
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values
t0=time.time()
TW = TW3()
Gsym = {g: col_Li3(sp.Rational(g), SS235) for g in Gargs}
sel = ['1/3', '2/3', '2/5', '-1/5', '1/10', '-5/3', '5', '3/5']
selq = [sp.Rational(g) for g in sel]
G8 = sp.Matrix.hstack(*[Gsym[g] for g in selq])
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
wlist = sorted(Gargs, key=lambda g:(nice(g), max(abs(Fr(g).numerator),abs(Fr(g).denominator)), str(g)))

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
li3vals = [Li3v(g) for g in selq]
pcolmap, pvalmap = {}, {}
def pcol(w,p):
    if (w,p) not in pcolmap: pcolmap[(w,p)] = col_prod(sp.Rational(w), p, SS235)
    return pcolmap[(w,p)]
def pval(w,p):
    if (w,p) not in pvalmap: pvalmap[(w,p)] = -Li2v(sp.Rational(w))*mp.log(p)
    return pvalmap[(w,p)]
z2 = mp.pi**2/6; z3 = mp.zeta(3)
def heights(PBspec):
    M = sp.Matrix.hstack(G8, *[pcol(w,p) for w,p in PBspec])
    if M.rank() < 27: return None
    out = []
    for e in (1,2,3):
        x = M.LUsolve(sp.Matrix(TW[e]))
        if M*x != sp.Matrix(TW[e]): return None
        c = [mp.mpf(str(sp.Rational(cc))) for cc in x]
        R = W3v[e] - mp.fsum([c[i]*li3vals[i] for i in range(8)]) - mp.fsum([c[8+j]*pval(*PBspec[j]) for j in range(19)])
        vec = [R, z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)]
        rel = mp.pslq(vec, tol=mp.mpf('1e-200'), maxcoeff=10**30, maxsteps=10000)
        if not rel: return None
        res = abs(mp.fsum([mp.mpf(r)*v for r,v in zip(rel,vec)]))
        out.append((R, rel, res))
    return out

PB_orig = [('-1/3',2),('-1/3',3),('-1/3',5),('-1/4',2),('-1/4',3),('-1/4',5),('1/4',2),('1/4',3),('1/4',5),
           ('-5/3',2),('-5/3',3),('-5/3',5),('-5/4',2),('-5/4',3),('-5/4',5),('1/5',3),('1/5',5),('5/8',2),('5/8',3)]
PB_orig = [(sp.Rational(w), p) for w,p in PB_orig]
h = heights(PB_orig)
print("ORIG PB:", [(mp.nstr(x[0],6), max(abs(r) for r in x[1]), mp.nstr(x[2],2)) for x in h], f"({time.time()-t0:.1f}s)")

random.seed(11)
best = (max(max(abs(r) for r in x[1]) for x in h), PB_orig)
for t in range(12):
    cur = []; curcols = sp.Matrix.hstack(G8, sp.zeros(27,0))
    tries = 0
    while len(cur) < 19 and tries < 800:
        tries += 1
        w = random.choice(wlist[:50]); p = random.choice((2,3,5))
        if (w,p) in cur: continue
        c = pcol(w,p)
        if sp.Matrix.hstack(curcols, c).rank() > curcols.shape[1]:
            cur.append((w,p)); curcols = sp.Matrix.hstack(curcols, c)
    if len(cur) < 19: continue
    hh = heights(cur)
    if hh is None: continue
    mx = max(max(abs(r) for r in x[1]) for x in hh)
    print(f"trial {t}: maxheight = {mx:.2e}  ({time.time()-t0:.1f}s)")
    if mx < best[0]: best = (mx, cur)
print("BEST height:", f"{best[0]:.2e}")
print("BEST PB:", [(str(w),p) for w,p in best[1]])
pickle.dump({'PB19': [(str(w),p) for w,p in best[1]]}, open("PB19_final.pkl","wb"))
