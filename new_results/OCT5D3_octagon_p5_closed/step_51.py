import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
import mpmath as mp
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import gengine2p

OCT5D = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(OCT5D+"parsed_mons.pkl","rb"))
rewrite = pickle.load(open("rewrite_all.pkl","rb"))
DPS = 400
G1, letters = gengine2p.gtab(430, 1440, cachedir="<project>/research/08_theory/OCT5D2_octagon_depth_test")
mp.mp.dps = DPS
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
Gv = {}
def gval(w):
    if w not in Gv: Gv[w] = G1[w2mpf(w)]
    return Gv[w]

# true weight parts: weight of term = |v| + e + f
Wt = {k: mp.mpf(0) for k in range(6)}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    cf = L2C**e * LU**f
    for n,(p,q) in cu.items():
        pf = mp.mpf(str(p)) + mp.mpf(str(q))*r2
        for w, gam in rewrite[(n,a,b,c,d)].items():
            k = len(w) + e + f
            if k <= 5:
                Wt[k] += pf*cf*q2f(gam)*gval(w)
Wt = {k: 960*v for k,v in Wt.items()}
for k in range(6):
    print("Wt%d =" % k, mp.nstr(Wt[k], 30))
print("sum:", mp.nstr(sum(Wt.values()), 30))
W = pickle.load(open("wparts_1650.pkl","rb"))
print("W5(deg5)  =", mp.nstr(mp.mpf(W['W5']), 30))
print("Wt5(top)  =", mp.nstr(Wt[5], 30))
print("difference:", mp.nstr(mp.mpf(W['W5'])-Wt[5], 10))
pickle.dump({f"Wt{k}": mp.nstr(Wt[k], 380) for k in range(6)}, open("wparts_true_430.pkl","wb"))
