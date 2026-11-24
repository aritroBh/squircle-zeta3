import sys; sys.path.insert(0,'.')
import sympy as sp, mpmath as mp
from d16_engine import build_B, Lambda, s_log, s_add, show
K=8
B=build_B(3,K)
L=s_add(Lambda(3,K), s_log(B,K))
show(L,'log((2-p)S/(4 sqrt3))')
mp.mp.dps=30
Z={k:mp.zeta(k) for k in range(2,K+1)}
for k in sorted(L):
    if k==0: continue
    val=mp.fsum(mp.mpf(sp.Rational(c).p)/sp.Rational(c).q*mp.fprod([Z[i] for i in m]) for m,c in L[k].items())
    print(f'tau^{k}: numeric {mp.nstr(val,15)}')
