"""a05_compare_B3.py -- exact comparison of the r-family HPL forms (a04, derived from the squircle integral)
with the weight-4 part of G1's fitted corner formula B3 (M3 = -z3 Phi + z2 Psi + R3).
Corner weight-4 part split by prefactor: F(c) = c^-1 F_-1 + F_0 + c F_1 (+ z2, z3 pieces from Psi, Phi).
At c = 1/sqrt2 the rational part is U^corner = F_0 and the sqrt2 part is V^corner = F_-1 + F_1/2."""
import json, pickle
from fractions import Fraction as Fr
from collections import defaultdict
G1 = "research/10_frontier/G1_zeta3_any_angle/"
R = pickle.load(open('a04_rfamily.pkl', 'rb'))
def parse(name):
    pre, w = name.split(' H')
    return int(pre.replace('c^', '')), tuple(int(x) for x in w.strip('[]').split(','))
corner = defaultdict(lambda: defaultdict(Fr))   # corner[e][(word, mono)]
for fname, mono in (("hpl_fit_k3_(0,0,0).json", ()), ("hpl_fit_k3_(1,0,0).json", ('z2',)),
                    ("hpl_fit_k3_(0,1,0).json", ('z3',))):
    sgn = -1 if mono == ('z3',) else 1          # M3 = ... - z3*Phi ; the json for z3 already IS the z3 coefficient
    for name, v in json.load(open(G1 + fname)):
        if name in ('1', 'c'):
            continue
        e, w = parse(name)
        # keep the weight-4 part only: |w| + weight(mono) = 4
        wt = len(w) + (2 if mono == ('z2',) else 3 if mono == ('z3',) else 0)
        if wt == 4:
            corner[e][(w, mono)] += Fr(v)
Uc = dict(corner[0]); Vc = defaultdict(Fr)
for k, v in corner[-1].items(): Vc[k] += v
for k, v in corner[1].items(): Vc[k] += v/2
print('prefactors in corner weight-4 part:', sorted(corner), ' c^-2 part:', dict(corner.get(-2, {})))
for name, A, B in (('U', R['U'], Uc), ('V', R['V'], dict(Vc))):
    keys = set(A) | set(B)
    diff = {k: A.get(k, 0) - B.get(k, 0) for k in keys if A.get(k, 0) != B.get(k, 0)}
    print(f'{name}: r-family terms {len(A)}, corner terms {len(B)}, exact difference terms: {len(diff)}')
    for k, v in sorted(diff.items(), key=str)[:40]:
        print('    ', k, v)
