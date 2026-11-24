import os, sys, time, pickle
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
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
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
ival = {k: sum(q2f(v)*gval(w) for w, v in r.items()) for k, r in rewrite.items()}

blocks = {'local': mp.mpf(0), 'single': mp.mpf(0), 'mixed': mp.mpf(0)}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    constfac = L2C**e * LU**f
    for n,(p,q) in cu.items():
        coef = (mp.mpf(str(p)) + mp.mpf(str(q))*r2)*constfac
        blocks[cls] += coef*ival[(n,a,b,c,d)]
blocks = {k: 960*v for k,v in blocks.items()}

for k in ['local','single','mixed']:
    tv = mp.mpf(targets[k])
    print(f"{k:7s} rewrite = {mp.nstr(blocks[k], 55)}")
    print(f"        target  = {mp.nstr(tv, 55)}")
    print(f"        residual = {mp.nstr(abs(blocks[k]-tv), 4)}")
tv = mp.mpf(targets['total'])
print("total   residual =", mp.nstr(abs(sum(blocks.values())-tv), 4))
quoted = mp.mpf("4.2299562389856442293498943256846704703130475386317")
print("mixed - quoted task value =", mp.nstr(blocks['mixed']-quoted, 4))
