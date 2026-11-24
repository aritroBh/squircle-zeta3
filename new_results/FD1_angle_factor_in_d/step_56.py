import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
from math import factorial
K=pickle.load(open(fd1+'/kern_sq_11_22.pkl','rb'))['kernels']
def q2f(t): return Fr(t[0],t[1])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
polemap = {((-1,1),(0,1),(0,1),(0,1)):'-a', ((1,1),(0,1),(0,1),(0,1)):'a',
           ((-3,1),(2,1),(0,1),(0,1)):'-b', ((0,1),(0,1),(0,1),(0,1)):'0'}
d=10
for lo in [11,12,13]:
    vals={}
    keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
    ok=True
    for key in keys:
        seqs = {}
        for k in range(lo, lo+10):
            for rho, cf in K[k][key]['real'].items():
                seqs.setdefault(rho, []).append(cf)
        for rho, s in seqs.items():
            if len(s) < d: continue
            p0=[q2f(cf[0]) for cf in s]
            t0 = fdiff(list(p0), d-1)
            if not all(x==t0[0] for x in t0): ok=False
            vals[(key,polemap.get(rho))]=t0[0]
    print(f"window k={lo}..{lo+9}: (10,0),-a = {vals.get(((10,0),'-a'))}, consistent={ok}")
