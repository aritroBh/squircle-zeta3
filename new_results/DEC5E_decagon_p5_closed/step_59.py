import sys, pickle
sys.path.insert(0,'.')
import mpmath as mp
from q5_lll import fit, report
mp.mp.dps = 800
phi = (1+mp.sqrt(5))/2
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
K3m = L3m + lf**3/6
li32 = pickle.load(open('li32_gold_800.pkl','rb'))
D32 = {}
D32['D32a'] = mp.mpf(li32['vals'][0]) + lf**5/120
D32['D32b'] = mp.mpf(li32['vals'][1]) + lf**5/15
D32['D32c'] = mp.mpf(li32['vals'][2]) + mp.mpf(23)/120*lf**5
D32['D32d'] = mp.mpf(li32['vals'][3]) + mp.mpf(73)/120*lf**5
for nm, v in D32.items():
    print(nm, "=", mp.nstr(v, 25))
INV = [('z2l2^3', z2*l2**3), ('z2l2^2lf', z2*l2**2*lf), ('z2l2lf^2', z2*l2*lf**2), ('z2lf^3', z2*lf**3),
       ('z3l2^2', z3*l2**2), ('z3l2lf', z3*l2*lf), ('z3lf^2', z3*lf**2), ('z2^2l2', z2**2*l2), ('z2^2lf', z2**2*lf),
       ('z5', z5), ('z2z3', z2*z3), ('z2L3m', z2*L3m), ('K5m', K5m), ('K5q', K5q),
       ('K4ml2', K4m*l2), ('K4mlf', K4m*lf), ('K4ql2', K4q*l2), ('K4qlf', K4q*lf),
       ('K3ml2^2', K3m*l2**2), ('K3ml2lf', K3m*l2*lf), ('K3mlf^2', K3m*lf**2)]
for nm, v in D32.items():
    c = fit(v, INV, cexp=int(800*0.90), nbest=2)
    print(f"--- {nm} in 21-atom span? top height %.2e resid %s" % (c[0][0], mp.nstr(mp.mpf(c[0][1]),2)))
