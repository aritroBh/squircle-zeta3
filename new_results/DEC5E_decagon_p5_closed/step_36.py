import sys, pickle
sys.path.insert(0,'../DEC5D_decagon_order5_depth')
import dec5d_lib as D
from fractions import Fraction as F
QS5 = D.QS5; ZERO, ONE = QS5(0), QS5(1)
P_LET = -D.PHI; Q_LET = D.PHI**2
rewrite = pickle.load(open('../DEC5D_decagon_order5_depth/rewrite_dec_exact.pkl','rb'))
parsed = pickle.load(open('../DEC5D_decagon_order5_depth/parsed_mons_dec.pkl','rb'))

# exact word values in basis {('1',a,b)=l2^a lf^b, ('z2',a,b)=z2*l2^a lf^b}
H = QS5(F(1,2)); T5 = QS5(F(2,5)); MT5 = QS5(F(-2,5)); F5 = QS5(F(1,5)); MF5 = QS5(F(-1,5))
VAL = {
 (): {('1',0,0): ONE},
 (ZERO,): {},
 (P_LET,): {('1',0,1): ONE},
 (Q_LET,): {('1',0,1): -ONE},
 (ZERO,ZERO): {},
 (ZERO,ONE): {('z2',0,0): -ONE},
 (ZERO,P_LET): {('z2',0,0): T5, ('1',0,2): -H},
 (ZERO,Q_LET): {('z2',0,0): MT5, ('1',0,2): ONE},
 (P_LET,ZERO): {('z2',0,0): MT5, ('1',0,2): H},
 (Q_LET,ZERO): {('z2',0,0): T5, ('1',0,2): -ONE},
 (P_LET,P_LET): {('1',0,2): H},
 (Q_LET,Q_LET): {('1',0,2): H},
 (P_LET,ONE): {('z2',0,0): MT5, ('1',0,2): ONE},
 (Q_LET,ONE): {('z2',0,0): T5, ('1',0,2): -H},
 (P_LET,Q_LET): {('z2',0,0): MF5, ('1',0,2): ONE},
 (Q_LET,P_LET): {('z2',0,0): F5, ('1',0,2): QS5(-2)},
}
def mul_ef(val, e, f):
    out = {}
    for (nm,a,b), c in val.items():
        out[(nm,a+e,b+f)] = out.get((nm,a+e,b+f), ZERO) + c
    return {k:v for k,v in out.items() if v != ZERO}
W = {0: {}, 1: {}, 2: {}}
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n,(p,q) in cu.items():
        cf = QS5(p, q)
        for w, gam in rewrite[(n,a,b,c,d)].items():
            k = len(w) + e + f
            if k > 2: continue
            for key, c2 in mul_ef(VAL[w], e, f).items():
                W[k][key] = W[k].get(key, ZERO) + cf*gam*c2
for k in range(3):
    print("Wt%d exact (x1200):" % k)
    for key, v in sorted(W[k].items()):
        print("   ", key, ":", 1200*v)
    print()
open('progress.log','a').write("step_33: EXACT Wt0=%r, Wt1=%r, Wt2=%r (QS5 arithmetic, proved word values)\n" % (
    {k: str(1200*v) for k,v in W[0].items()}, {k: str(1200*v) for k,v in W[1].items()}, {k: str(1200*v) for k,v in W[2].items()}))
