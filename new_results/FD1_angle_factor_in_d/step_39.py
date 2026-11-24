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
data = pickle.load(open(f'{fd1}/kern20_c12.pkl','rb'))
K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
d=10
keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
for key in keys:
    seqs = {}
    for k in range(d, 21):
        P = K[k][key]
        for rho, cf in P['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
    for rho, s in seqs.items():
        top = fdiff(list(s[1:]), d-1)
        beta = top[0]/factorial(d-1)
        pole = Fr(*rho)
        polename = '-b' if pole==-r else str(int(pole))
        sb = squircle_beta(d, key[0], key[1], {-1:'-a',1:'a',0:'0'}.get(int(pole),'-b') if pole!=-r else '-b')
        print(key, polename, "beta=", beta, " sq=", sb)
