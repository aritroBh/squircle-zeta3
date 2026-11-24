"""step_22.py -- weight-three slope kernels: residues of dK/dk (constant for k>=3), at 4 angles."""
import sys, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr

def qpart(t):
    return Fr(*t[0]), Fr(*t[1])

ANGLES = ['c12', 'c35', 'c13', 'c45']
for name in ANGLES:
    data = pickle.load(open(f'kern_{name}.pkl', 'rb'))
    K = data['kernels']
    r = Fr(*data['r'])
    print(f"=== {name} (c={Fr(*data['c'])}, r={r}) ===")
    for key in ((2, 0), (1, 1)):
        # slope = K[k+1] - K[k] for k = 4..10 (should be constant)
        slopes = {}
        for k in range(4, 11):
            P1, P0 = K[k+1][key], K[k][key]
            poles = set(P1['real']) | set(P0['real'])
            for rho in poles:
                c1 = qpart(P1['real'].get(rho, ((0,1),(0,1))))[0] if rho in P1['real'] else Fr(0)
                c0 = qpart(P0['real'].get(rho, ((0,1),(0,1))))[0] if rho in P0['real'] else Fr(0)
                slopes.setdefault(rho, []).append(c1 - c0)
        for rho, vals in sorted(slopes.items(), key=str):
            const = all(v == vals[0] for v in vals)
            print(f"  key{key} pole{rho}: slope = {vals[0]}  (constant over k=4..10: {const})")
print("STEP22 DONE")
