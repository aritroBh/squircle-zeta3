"""step_31.py -- top Newton table beta_d(r; c): leading coefficients of residue polynomials in k.
beta_d^{(ij)}(rho) = [k^{d-1}] res_rho(K^{(k)}_{ij}) = (d-1)-th finite difference / (d-1)!.
Extract exactly at c=1/2, 3/5, 1/3, 4/5 (and note the pattern)."""
import sys, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from math import factorial

def qpart(t):
    return Fr(*t[0]), Fr(*t[1])

def fdiff(seq, n):
    for _ in range(n):
        seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq

ANGLES = ['c12', 'c35', 'c13', 'c45']
tables = {}
for name in ANGLES:
    data = pickle.load(open(f'kern_{name}.pkl', 'rb'))
    K = data['kernels']
    r = Fr(*data['r'])
    tables[name] = {}
    for d in range(1, 7):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, 15):
                P = K[k][key]
                for rho, cf in P['real'].items():
                    seqs.setdefault(rho, []).append(qpart(cf)[0])
            for rho, s in seqs.items():
                # s[k0] = res at k=d+k0; poly for k>=d+1 -> use s[1:]; leading coeff = Delta^{d-1} s[1:]/(d-1)!
                top = fdiff(list(s[1:]), d-1)
                beta = top[0]/factorial(d-1)
                if beta != 0:
                    tables[name][(d, key, rho)] = beta
    print(f"=== {name} (r={r}) ===")
    for d in range(1, 7):
        entries = {kk: vv for kk, vv in tables[name].items() if kk[0] == d}
        if not entries: continue
        print(f"  d={d}:")
        for (dd, key, rho), beta in sorted(entries.items()):
            print(f"     key{key} pole{rho}: beta = {beta}")
import pickle as pk
pk.dump(tables, open('beta_tables.pkl', 'wb'))
print("STEP31 DONE")
