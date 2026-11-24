import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.getcwd())
from chk_feed import *
from fastred import pf_to_sympy, v
import sympy as sp
from fractions import Fraction as Fr
def lagr_vec(xs, ys, x):
    s=[Fr(0)]*4
    for i,xi in enumerate(xs):
        L=Fr(1)
        for j,xj in enumerate(xs):
            if i!=j: L*=Fr(x-xj, xi-xj)
        for t in range(4): s[t]+=ys[i][t]*L
    return tuple(s)
Dl={}
for dl in range(9):
    P=dl//2; xs=list(range(7,7+P+1)); ys=[H[(3,3)][dl][k] for k in xs]
    st=lagr_vec(xs,ys,6)
    Dl[dl]=tuple(a-b for a,b in zip(H[(3,3)][dl][6],st))
def Psi(dls):
    hD={10-dl:(QSI(*Dl[dl]) if dl in dls else QSI(0)) for dl in range(9)}
    G=gm_from_tower(hD)
    expr=pf_to_sympy(G)
    return sp.simplify(sp.together(expr-expr.subs(v,1/v)))
full=Psi(set(range(9)))
print('Psi_Delta(v)=Gm(v)-Gm(1/v) for the full discrepancy simplifies to:', full)
for dl in range(9):
    r=Psi({dl})
    print('single depth',dl,'Psi == 0 ?', r==0)
