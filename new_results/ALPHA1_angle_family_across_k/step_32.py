"""step_32.py -- squircle beta table from published kern_k.pkl (k<=18), then ratios master_d(c)/master_d(squircle)."""
import sys, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from math import factorial

TH = '<project>/research/08_theory/WTS_weight_k_structure/theory/'
def qsi(t):
    return tuple(Fr(x) for x in t)

def fdiff(seq, n):
    for _ in range(n):
        seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq

# squircle kernels k<=18
K = {}
import os
for k in range(1, 19):
    fn = TH + f'kern_{k}.pkl'
    if os.path.exists(fn):
        K[k] = pickle.load(open(fn, 'rb'))
print("squircle levels available:", sorted(K))

def beta_table(K, pole_map, dmax=8, kmax=None):
    """pole_map: dict from my pole label ('0','1','-1','r','-r') to squircle pkl pole key."""
    tab = {}
    if kmax is None: kmax = max(K)
    for d in range(1, dmax+1):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, kmax+1):
                if k not in K or key not in K[k]: continue
                P = K[k][key]
                for rho, cf in P['real'].items():
                    seqs.setdefault(rho, []).append(qsi(cf))
            for rho, s in seqs.items():
                if len(s) < d+2: continue
                for comp in range(4):
                    sc = [x[comp] for x in s]
                    top = fdiff(list(sc[1:]), d-1)
                    beta = top[0]/factorial(d-1)
                    if beta != 0:
                        tab.setdefault((d, key, rho), [Fr(0)]*4)[comp] = beta
    return tab

# squircle pkl poles: (p0,p1) = p0+p1 sqrt2 ; a=(1,1), b=(-1,1), 0=(0,1)
A2 = (Fr(1), Fr(1)); B2 = (Fr(-1), Fr(1)); Z2 = (Fr(0), Fr(0))
sqt = beta_table(K, None, dmax=8, kmax=18)
# map to w labels
def wlab(rho):
    if rho == A2: return '1'
    if rho == tuple(-x for x in A2): return '-1'
    if rho == B2: return 'r'
    if rho == tuple(-x for x in B2): return '-r'
    if rho == Z2: return '0'
    return str(rho)
sqtab = {}
for (d, key, rho), beta in sqt.items():
    # keep only rational part (p2=p3=0 expected)
    beta = sqt[(d, key, rho)]
    assert beta[2] == 0 and beta[3] == 0, (d, key, rho, beta)
    sqtab[(d, key, wlab(rho))] = (beta[0], beta[1])
print("squircle beta table (rational, sqrt2 parts):")
for d in range(1, 9):
    entries = {kk: vv for kk, vv in sqtab.items() if kk[0] == d}
    if not entries: continue
    print(f"  d={d}:")
    for (dd, key, rho), beta in sorted(entries.items()):
        print(f"     key{key} pole {rho}: {beta}")
pickle.dump(sqtab, open('beta_squircle.pkl', 'wb'))
print("STEP32 DONE")
