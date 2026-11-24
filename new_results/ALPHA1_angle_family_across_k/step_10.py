"""step_10.py -- degree analysis (exact degree d-1? k=d on the polynomial?), weight-two residues."""
import sys, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr

def qpart(t):
    return Fr(*t[0]), Fr(*t[1])

def fdiff(seq, n):
    for _ in range(n):
        seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq

ANGLES = ['c12', 'c35', 'c13', 'c45']
print("== exact degree check: is Delta^{d-1} res nonzero for some pole of each key? and is k=d ON the poly? ==")
for name in ANGLES:
    data = pickle.load(open(f'kern_{name}.pkl', 'rb'))
    K = data['kernels']
    print(f"=== {name} ===")
    for d in range(1, 6):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, 15):
                P = K[k][key]
                for rho, cf in P['real'].items():
                    seqs.setdefault(rho, []).append(qpart(cf)[0])
            maxdeg = -1; kd_on = True
            for rho, s in seqs.items():
                # s[k0] = res at k = d+k0
                # degree of polynomial for k>=d+1: largest n with Delta^n (s[1:]) != 0
                for n in range(0, d):
                    if any(x != 0 for x in fdiff(list(s[1:]), n)) : maxdeg = n
                # is k=d on the polynomial? extrapolate poly from k=d+1.. back to k=d:
                # value at k=d predicted by poly = (-1)^d * (Delta^d s with s[d] free)... use Newton:
                sp = s[1:]   # k=d+1..14
                # predicted s[0] via finite difference: s[0] = sum_{j} (-1)^{j+1} ... use Lagrange at first d points
                pts = sp[:d]  # k = d+1 .. 2d
                # Lagrange extrapolation to k=d
                pred = Fr(0)
                for j in range(d):
                    term = pts[j]
                    for m in range(d):
                        if m != j:
                            term *= Fr(d - (d+1+m), (d+1+j) - (d+1+m))
                    pred += term
                if pred != s[0]:
                    kd_on = False
            print(f"  key{key}: max res degree for k>=d+1 = {maxdeg} (expect {d-1}); k=d ON poly: {kd_on}")
print()
print("== WEIGHT TWO: K^{(k)}_{10} residues for k>=2 (should be k-independent) ==")
for name in ANGLES:
    data = pickle.load(open(f'kern_{name}.pkl', 'rb'))
    K = data['kernels']
    r = Fr(*data['r'])
    print(f"=== {name} (r={r}) ===")
    seen = {}
    for k in (2, 3, 5, 8):
        P = K[k][(1, 0)]
        res = {rho: qpart(cf)[0] for rho, cf in P['real'].items()}
        seen[k] = res
    poles = sorted(seen[2], key=str)
    for rho in poles:
        vals = [seen[k].get(rho) for k in (2,3,5,8)]
        rn, rd = rho
        print(f"  pole ({rn}/{rd}): residues k=2,3,5,8 = {vals}")
print("STEP10 DONE")
