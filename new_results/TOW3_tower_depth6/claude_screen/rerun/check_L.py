import sys, time; sys.path.insert(0,'.')
import sympy as sp
import tow2lib as T, tow3lib as T3
d, m = T.d, T.m
Lp = T3._load_L()
print('L_parsed keys:', sorted(Lp.keys()))
print('free symbols of L_parsed[r]:', [sorted(map(str, Lp[r].free_symbols)) for r in sorted(Lp)][:3])
pts = [(sp.Rational(7), sp.Rational(3)), (sp.Rational(11,3), sp.Rational(5,2)), (sp.Rational(23), sp.Rational(17))]
t0=time.time()
for r in range(1, 7):
    res = []
    for (dv, mv) in pts:
        a = T.L(r).subs({d: dv, m: mv}); b = Lp[r].subs({d: dv, m: mv})
        diff = sp.simplify(sp.expand(a - b))
        res.append(diff)
    print(f'r={r}: L(r) [tow2lib proved assembly from N_1..N_4] - L_parsed[r] at 3 exact points:', res, round(time.time()-t0,1), flush=True)
# a5, a6 via the tow2lib proved Bell assembly vs tow3lib.arho at points
for dl in (5, 6):
    res = []
    for (dv, mv) in pts:
        a = T.complete_bell(dl).subs({d: dv, m: mv}); b = T3.arho(dl).subs({d: dv, m: mv})
        res.append(sp.simplify(sp.expand(a - b)))
    print(f'a_{dl}: tow2lib.complete_bell - tow3lib.arho at 3 exact points:', res, round(time.time()-t0,1), flush=True)
