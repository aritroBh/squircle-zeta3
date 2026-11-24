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
polemap = {((1,1),(-1,1),(0,1),(0,1)):'-a', ((-1,1),(1,1),(0,1),(0,1)):'a',
           ((7,1),(-5,1),(0,1),(0,1)):'-b', ((0,1),(0,1),(0,1),(0,1)):'0'}
for d in [8,9,10]:
    keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
    print(f"--- d={d} ---")
    for key in keys:
        seqs = {}
        for k in range(d, 21):
            if k not in K or key not in K[k]: continue
            for rho, cf in K[k][key]['real'].items():
                seqs.setdefault(rho, []).append(cf)
        for rho, s in seqs.items():
            if len(s) < d+1: continue
            vals=[q2f(cf[0]) for cf in s]  # p0 part; assert real
            assert all(q2f(cf[1])==0 and q2f(cf[2])==0 and q2f(cf[3])==0 for cf in s)
            top = fdiff(list(vals[1:]), d-1)
            assert all(x==top[0] for x in top), (d,key,rho)
            Dtop = top[0]*factorial(d-1)/factorial(d-1)  # beta*(d-1)! = D
            beta = top[0]/factorial(d-1)
            D = beta*factorial(d-1)
            pole = polemap.get(rho, ('?',rho))
            cf_form = squircle_D(d, key[0], key[1], pole) if pole in ('-a','a','-b','0') else None
            flag = "" if cf_form==D else f"  <-- closed form {cf_form} ratio {D/cf_form if cf_form else 'NA'}"
            if D!=0 or (cf_form not in (None,Fr(0))): print(f"  {key} {pole}: D={D}{flag}")
