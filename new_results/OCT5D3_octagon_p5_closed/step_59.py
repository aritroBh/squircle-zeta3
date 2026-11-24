import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import mpmath as mp
OCT5D = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(OCT5D+"parsed_mons.pkl","rb"))
rewrite = pickle.load(open("rewrite_all.pkl","rb"))

# G_*((a,);1) values as exact QS2 multiples of lam: (0)->0, (1)->0, (p)->-1/2, (q)->+1/2
GLAM = {ZERO: Q(0), ONE: Q(0), P_LET: Q(-1,0)/2, Q_LET: Q(1,0)/2}
W1_lam = Q(0); W1_L = Q(0)
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n,(p,q) in cu.items():
        pf = Q(p,q)
        for w, gam in rewrite[(n,a,b,c,d)].items():
            k = len(w)+e+f
            if k != 1: continue
            if len(w) == 0:   # then e+f = 1: constant log
                if e == 1: W1_lam += pf*gam
                else: W1_L += pf*gam
            else:             # len(w)=1, e=f=0
                W1_lam += pf*gam*GLAM[w[0]]
W1_lam *= 960; W1_L *= 960
print("Wt1 = (%s)*lam + (%s)*L" % (W1_lam, W1_L))
mp.mp.dps = 100
r2 = mp.sqrt(2)
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*r2
val = q2f(W1_lam)*mp.log(2) + q2f(W1_L)*mp.log(1+r2)
Wt1 = mp.mpf(pickle.load(open("wparts_true_1650.pkl","rb"))['Wt1'])
print("exact-form value:", mp.nstr(val, 40))
print("numeric Wt1     :", mp.nstr(Wt1, 40))
print("difference:", mp.nstr(abs(val-Wt1), 3))
