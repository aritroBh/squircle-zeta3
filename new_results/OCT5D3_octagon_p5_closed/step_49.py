import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
import mpmath as mp
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import gengine2p

t0 = time.time()
OCT5D = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(OCT5D+"parsed_mons.pkl","rb"))
# exact rewrites for ALL monomials (degrees 0..5)
need = set()
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n in cu: need.add((n,a,b,c,d))
if os.path.exists("rewrite_all.pkl"):
    rewrite = pickle.load(open("rewrite_all.pkl","rb"))
else:
    rewrite = {k: E.monomial_integral(*k) for k in sorted(need)}
    pickle.dump(rewrite, open("rewrite_all.pkl","wb"))
print("rewrites:", len(rewrite), "(%d integrals)" % len(need), "%.1fs" % (time.time()-t0))

DPS = 1650
G1, letters = gengine2p.gtab(DPS, 5700, cachedir=".")
mp.mp.dps = DPS
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
Gv = {}
def gval(w):
    if w not in Gv: Gv[w] = G1[w2mpf(w)]
    return Gv[w]

ival = {k: sum(q2f(v)*gval(w) for w, v in r.items()) for k, r in rewrite.items()}
W = {k: mp.mpf(0) for k in range(6)}
for mon, cu in parsed:
    k = sum(mon)
    a,b,c,d,e,f = mon
    cf = L2C**e * LU**f
    for n,(p,q) in cu.items():
        W[k] += (mp.mpf(str(p)) + mp.mpf(str(q))*r2)*cf*ival[(n,a,b,c,d)]
W = {k: 960*v for k,v in W.items()}
for k in range(6):
    print("W%d =" % k, mp.nstr(W[k], 40))
tot = sum(W.values())
print("sum  =", mp.nstr(tot, 40))
print("assembly %.1fs" % (time.time()-t0))
pickle.dump({f"W{k}": mp.nstr(W[k], DPS-20) for k in range(6)}, open("wparts_1650.pkl","wb"))

# cross-checks
p5 = pickle.load(open(OCT5D+"p5_800.pkl","rb"))
P5over_u_800 = mp.mpf(p5['P5over_u'])
print("P5/u (800) - sum(Wk):", mp.nstr(abs(P5over_u_800-tot), 3))
w5_430 = mp.mpf(pickle.load(open("<project>/research/08_theory/OCT5D2_octagon_depth_test/w5_430.pkl","rb"))['W5'])
print("W5 - w5_430:", mp.nstr(abs(W[5]-w5_430), 3))
