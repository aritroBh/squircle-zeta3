import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 250
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import gengine

SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(SRC+"parsed_mons.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
G1, letters = gengine.gtab(250, 830)
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
Gcache = {}
def gval(w):
    if w not in Gcache: Gcache[w] = G1[w2mpf(w)]
    return Gcache[w]
t0=time.time()
ival = {k: sum(q2f(v)*gval(w) for w, v in r.items()) for k, r in rewrite.items()}
blocks = {'local': mp.mpf(0), 'single': mp.mpf(0), 'mixed': mp.mpf(0)}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    cf = L2C**e * LU**f
    for n,(p,q) in cu.items():
        blocks[cls] += (mp.mpf(str(p)) + mp.mpf(str(q))*r2)*cf*ival[(n,a,b,c,d)]
blocks = {k: 960*v for k,v in blocks.items()}
W5 = sum(blocks.values())
print("assembled in %.1fs" % (time.time()-t0))
print("W5/u (250 digits) =", mp.nstr(W5, 60))
targets = pickle.load(open("targets110.pkl","rb"))
print("vs 110-digit target: ", mp.nstr(abs(W5-mp.mpf(targets['total'])), 4))
print("mixed block:         ", mp.nstr(abs(blocks['mixed']-mp.mpf(targets['mixed'])), 4))
pickle.dump({'W5': mp.nstr(W5, 240), 'blocks': {k: mp.nstr(v,240) for k,v in blocks.items()}},
            open("w5_250.pkl","wb"))
