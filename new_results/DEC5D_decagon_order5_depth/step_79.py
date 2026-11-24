
import pickle
from fractions import Fraction as F
import collections
import dec5d_lib as D
import mpmath as mp
mp.mp.dps = 250
QS5 = D.QS5
exactI = {k: {n: F(v) for n, v in vv.items()} for k, vv in
          ((k, pickle.load(open("localI_exact.pkl","rb"))[k]) for k in pickle.load(open("localI_exact.pkl","rb")))}
parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
# assemble: local block = 1200 * sum_{(a,b,0,0,e,f), sum=5} (p+q s5) l2^e lf^f I
# accumulate exact coefficients for basis atoms: atom = (zeta-name, e, f) with coeff in QS5
acc = collections.defaultdict(lambda: QS5(0))
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    if not (c==0 and d==0) or sum(mon)!=5: continue
    for n,(p,q) in cu.items():
        cf = QS5(F(p), F(q))
        for zn, zc in exactI[(n,a,b)].items():
            acc[(zn, e, f)] += cf*zc
acc = {k: 1200*v for k, v in acc.items() if v != QS5(0)}
print("exact local block decomposition (atom, l2^e lf^f) : coefficient in Q(sqrt5)")
for k in sorted(acc):
    print("   %-6s l2^%d lf^%d : %s" % (k[0], k[1], k[2], acc[k]))
# numeric check
r5 = mp.sqrt(5); phi = (1+r5)/2
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
ZV = {'1': mp.mpf(1), 'z2': z2, 'z3': z3, 'z2^2': z2**2, 'z5': z5, 'z2z3': z2*z3}
val = sum(float(acc[k])*ZV[k[0]]*l2**k[1]*lf**k[2] for k in acc)
local250 = mp.mpf(pickle.load(open("local250.pkl","rb")))
print("\nnumeric local block:", mp.nstr(val, 55))
print("quad local block :", mp.nstr(local250, 55))
print("residual:", mp.nstr(abs(val-local250), 3))
# zeta(5) and zeta(2)zeta(3) coefficients
print("\nzeta(5) coeff:", acc.get(('z5',0,0)), " paper: (73125/64-30825/64*sqrt5)")
print("zeta(2)zeta(3) coeff:", acc.get(('z2z3',0,0)), " paper: (-64875/64+32775/64*sqrt5)")
pickle.dump(acc, open("local_block_exact.pkl","wb"))
