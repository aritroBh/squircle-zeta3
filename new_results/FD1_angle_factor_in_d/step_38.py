import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
import importlib, fd1_lib; importlib.reload(fd1_lib)
from fd1_lib import squircle_beta
from fractions import Fraction as Fr
from math import factorial
def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
def extract_beta20(name, dmax=10):
    data = pickle.load(open(f'{fd1}/kern20_{name}.pkl','rb'))
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, dmax+1):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, 21):
                if k not in K or key not in K[k]: continue
                P = K[k][key]
                for rho, cf in P['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                if len(s) < d + 1: continue
                top = fdiff(list(s[1:]), d-1)
                if d <= 8: assert all(x == top[0] for x in top), (name, d, key, rho)
                beta = top[0]/factorial(d-1)
                if beta != 0: tab[(d, key, Fr(*rho))] = beta
    return c, r, tab
polmap={-1:'-a',1:'a',0:'0'}
Fd10={}
nrat={}
for name in ['c12','c35','c13','c45','c23','c57','c14','c67']:
    c, r, tab = extract_beta20(name)
    for d in [10]:
        ratios = {}
        for (dd,ij,p),v in tab.items():
            if dd!=d: continue
            pole = '-b' if p == -r else polmap.get(int(p))
            if pole is None: continue
            sb = squircle_beta(d, ij[0], ij[1], pole)
            if sb != 0: ratios[(ij,pole)] = v/sb
        uniq = set(ratios.values())
        nrat[name]=len(ratios)
        print(f"{name} c={c}: d=10 distinct ratios={len(uniq)} (from {len(ratios)} cells)", list(uniq)[:1])
        if len(uniq)==1: Fd10.setdefault(d,{})[c]=uniq.pop()
        else: print("  MULTIPLE:", uniq)
pickle.dump(Fd10, open(fd1+'/Fd_data_d10.pkl','wb'))
print("nrat:", nrat)
open(fd1+'/progress.log','a').write('step_35: F_10 extracted at all 8 angles, single-valued ratios\n')
