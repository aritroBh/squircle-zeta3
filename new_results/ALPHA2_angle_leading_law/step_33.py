import pickle
from fractions import Fraction as Fr
from math import factorial

def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n):
        seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq

def extract_beta(name, kmax):
    data = pickle.load(open(f'<project>/research/08_theory/ALPHA2_angle_leading_law/kern_{name}.pkl','rb'))
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, 8):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, kmax+1):
                P = K[k][key]
                for rho, cf in P['real'].items():
                    seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                # need k>=d+1 polynomial: use s[1:], need at least d+2 points for 1 spare
                if len(s) < d + 3: continue
                top = fdiff(list(s[1:]), d-1)
                assert all(x == top[0] for x in top), (name, d, key, rho, top[:4])
                beta = top[0]/factorial(d-1)
                if beta != 0:
                    tab[(d, key, rho)] = beta
    return c, r, tab

for name in ['c23','c57','c14']:
    c, r, tab = extract_beta(name, 16)
    print(f"== {name}: c={c}, r={r}, entries={len(tab)}")
pickle.dump('ok', open('/tmp/ok.pkl','wb'))
print("extraction works")
