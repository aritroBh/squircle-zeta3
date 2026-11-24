import sys, time
sys.dont_write_bytecode=True
sys.path.insert(0,'<project>/research/08_theory/TOW4_tower_depth8')
from core import *
import sympy as sp
import tow2lib as T2
# (c) R_l from definition vs N_l/Pi_l
d_, m_ = T2.d, T2.m
bad=0; cells=0; t0=time.time()
for d in range(1, 41):
    for m in range(0, d+1):
        c = S_taylor(d, m, 4)
        for l in range(1,5):
            Rl = c[l]*factorial(l)/c[0]          # S^(l)(t*)/S(t*)
            exg = sp.expand(T2.R[l].subs({d_: sp.Integer(d), m_: sp.Integer(m)}))
            re_, im_ = sp.re(exg), sp.im(exg)
            re_=sp.Rational(re_); im_=sp.Rational(im_)
            cells+=1
            if not (Rl.r == fmpq(int(re_.p),int(re_.q)) and Rl.i == fmpq(int(im_.p),int(im_.q))): bad+=1
print(f'R_l (l=1..4) from the S_m definition vs N_l/Pi_l: all 1<=d<=40, 0<=m<=d: cells={cells} bad={bad} ({time.time()-t0:.0f}s)')
