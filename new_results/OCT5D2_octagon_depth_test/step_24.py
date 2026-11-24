import os, sys, time, pickle, math
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 170
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import gengine

SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(SRC+"parsed_mons.pkl","rb"))
targets = pickle.load(open("targets110.pkl","rb"))

t0 = time.time()
# collect needed (n,a,b,c,d) per class at weight 5
need = {'local': set(), 'single': set(), 'mixed': set()}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    for n in cu: need[cls].add((n,a,b,c,d))
print({k: len(v) for k,v in need.items()}, "total", sum(map(len, need.values())))

# exact rewrite of each needed integral (cache)
rewrite = {}
for cls in need:
    for (n,a,b,c,d) in sorted(need[cls]):
        rewrite[(n,a,b,c,d)] = E.monomial_integral(n,a,b,c,d)
print("exact rewrite done in %.1fs, %d integrals" % (time.time()-t0, len(rewrite)))
nwords = set()
for r in rewrite.values(): nwords.update(r.keys())
print("distinct G-words appearing:", len(nwords), "max weight:", max(map(len, nwords)))
# check: no word starting with 1 has nonzero coefficient
bad1 = [(k,w) for k,r in rewrite.items() for w,v in r.items() if w and w[0]==ONE and v!=ZERO]
print("words starting with 1 with nonzero coeff:", len(bad1))
# save exact rewrite (QS2 pairs of Fractions are picklable)
pickle.dump(rewrite, open("rewrite_exact.pkl","wb"))

# numeric evaluation against G-table
G1, letters = gengine.gtab(170, 480)
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
L0, L1, LP, LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w):
    out = []
    for x in w:
        if x == ZERO: out.append(L0)
        elif x == ONE: out.append(L1)
        elif x == P_LET: out.append(LP)
        elif x == Q_LET: out.append(LQ)
        else: raise ValueError(x)
    return tuple(out)
Gcache = {}
def gval(w):
    if w not in Gcache: Gcache[w] = G1[w2mpf(w)]
    return Gcache[w]

ival = {}
for k, r in rewrite.items():
    ival[k] = sum(q2f(v)*gval(w) for w, v in r.items())

# spot-check individual integrals against itab110
ITAB = {k: mp.mpf(v) for k, v in pickle.load(open(SRC+"itab110.pkl","rb")).items()}
worst = mp.mpf(0); worstk=None
for k, v in ival.items():
    n,a,b,c,d = k
    ref = ITAB[(n,a,b,c,d)]
    if abs(v-ref) > worst: worst, worstk = abs(v-ref), k
print("worst single-integral residual vs itab110:", mp.nstr(worst, 4), "at", worstk)
