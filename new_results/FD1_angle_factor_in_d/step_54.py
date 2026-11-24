import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
import importlib, fd1_lib; importlib.reload(fd1_lib)
from fd1_lib import squircle_D
from fractions import Fraction as Fr
from math import factorial
K=pickle.load(open(fd1+'/kern_sq_11_20.pkl','rb'))['kernels']
def q2f(t): return Fr(t[0],t[1])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
polemap = {((-1,1),(0,1),(0,1),(0,1)):'-a', ((1,1),(0,1),(0,1),(0,1)):'a',
           ((-3,1),(2,1),(0,1),(0,1)):'-b', ((0,1),(0,1),(0,1),(0,1)):'0'}
d=10
keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
print(f"--- squircle d={d} top (TRUE) ---")
for key in keys:
    seqs = {}
    for k in range(11, 21):
        for rho, cf in K[k][key]['real'].items():
            seqs.setdefault(rho, []).append(cf)
    for rho, s in seqs.items():
        if len(s) < d: continue
        p0=[q2f(cf[0]) for cf in s]; p1=[q2f(cf[1]) for cf in s]
        t0 = fdiff(list(p0), d-1); t1 = fdiff(list(p1), d-1)
        assert all(x==t0[0] for x in t0) and all(x==t1[0] for x in t1), (key,rho)
        D0, D1 = t0[0], t1[0]
        pole = polemap.get(rho)
        cf_form = squircle_D(d, key[0], key[1], pole)
        flag = "" if (cf_form==D0 and D1==0) else f"  <== closed {cf_form}, ratio={D0/cf_form if cf_form else 'NA'}, sq2={D1}"
        if D0!=0 or D1!=0: print(f"  {key} {pole}: D={D0}{flag}")
