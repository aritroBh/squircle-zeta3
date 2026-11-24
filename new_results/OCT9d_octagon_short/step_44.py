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

DPS = 80
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
pcolmap = {}
pvalmap = {}
def pcol(w,p):
    if (w,p) not in pcolmap: pcolmap[(w,p)] = col_prod(sp.Rational(w), p, SS235)
    return pcolmap[(w,p)]
def pval(w,p):
    if (w,p) not in pvalmap: pvalmap[(w,p)] = -Li2v(sp.Rational(w))*mp.log(p)
    return pvalmap[(w,p)]

def evalPB(PBspec, dps=80):
    mp.mp.dps = dps
    M = sp.Matrix.hstack(G8, *[pcol(w,p) for w,p in PBspec])
    if M.rank() < 27: return None
    out = []
    for e in (1,2,3):
        x = M.LUsolve(sp.Matrix(TW[e]))
        if M*x != sp.Matrix(TW[e]): return None
        c = [mp.mpf(str(sp.Rational(cc))) for cc in x]
        R = W3v[e] - mp.fsum([c[i]*li3vals[i] for i in range(8)]) - mp.fsum([c[8+j]*pval(*PBspec[j]) for j in range(19)])
        out.append(R)
    return out

def height_of_R(R, dps=200):
    mp.mp.dps = dps
    z2 = mp.pi**2/6; z3 = mp.zeta(3)
    vec = [R, z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)]
    rel = mp.pslq(vec, tol=mp.mpf('1e-170'), maxcoeff=10**30, maxsteps=10000)
    if not rel: return None, None
    res = abs(mp.fsum([mp.mpf(r)*v for r,v in zip(rel,vec)]))
    return max(abs(r) for r in rel), res

random.seed(5)
def random_PB():
    cur = []; curcols = sp.Matrix.hstack(G8, sp.zeros(27,0))
    tries = 0
    while len(cur) < 19 and tries < 500:
        tries += 1
        w = random.choice(wlist[:60]); p = random.choice((2,3,5))
        if (w,p) in cur: continue
        c = pcol(w,p)
        if sp.Matrix.hstack(curcols, c).rank() > curcols.shape[1]:
            cur.append((w,p)); curcols = sp.Matrix.hstack(curcols, c)
    return cur if len(cur)==19 else None

best = None
for t in range(14):
    PB = random_PB()
    if PB is None: continue
    Rs = evalPB(PB)
    if Rs is None: continue
    score = max(abs(R) for R in Rs)
    hh = height_of_R(Rs[0])  # e1 constant height as second check
    print(f"trial {t}: max|R| = {mp.nstr(score, 3)}  e1 height = {hh[0]}  ({time.time()-t0:.1f}s)")
    if best is None or score < best[0]:
        best = (score, PB)
print("BEST max|R|:", mp.nstr(best[0],3))
print("BEST PB:", [(str(w),p) for w,p in best[1]])
pickle.dump({'PB19_rand': [(str(w),p) for w,p in best[1]]}, open("PB19_rand.pkl","wb"))
