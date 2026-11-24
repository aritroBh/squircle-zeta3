import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
import mpmath as mp
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
t0 = time.time()
OCT5D = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(OCT5D+"parsed_mons.pkl","rb"))
rewrite = pickle.load(open("rewrite_all.pkl","rb"))
raw, letters = pickle.load(open("gtab2_d1650_N5700.pkl","rb"))
mp.mp.dps = 1650
G1 = {k: mp.mpf(v) for k, v in raw.items()}
print("loaded %.1fs" % (time.time()-t0))
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1+r2)
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
def w2mpf(w): return tuple(L0 if x==ZERO else L1 if x==ONE else LP if x==P_LET else LQ for x in w)
Gv = {}
def gval(w):
    if w not in Gv: Gv[w] = G1[w2mpf(w)]
    return Gv[w]
# exact Wt0 check (QS2 arithmetic) + numeric parts
Wt0_exact = Q(0)
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n,(p,q) in cu.items():
        pf = Q(p, q)
        for w, gam in rewrite[(n,a,b,c,d)].items():
            if len(w)+e+f == 0:
                Wt0_exact += pf*gam
print("Wt0 exact =", Wt0_exact*960)
t0 = time.time()
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
print("assembled %.1fs" % (time.time()-t0))
for k in range(6):
    print("Wt%d =" % k, mp.nstr(Wt[k], 35))
print("sum:", mp.nstr(sum(Wt.values()), 35))
p5 = pickle.load(open(OCT5D+"p5_800.pkl","rb"))
print("P5/u(800) - sum:", mp.nstr(abs(mp.mpf(p5['P5over_u'])-sum(Wt.values())), 3))
w5_430 = pickle.load(open("wparts_true_430.pkl","rb"))
for k in range(1,6):
    print(f"Wt{k} vs 430-digit:", mp.nstr(abs(Wt[k]-mp.mpf(w5_430[f'Wt{k}'])), 3))
pickle.dump({f"Wt{k}": mp.nstr(Wt[k], 1620) for k in range(6)}, open("wparts_true_1650.pkl","wb"))
print("saved wparts_true_1650.pkl")
