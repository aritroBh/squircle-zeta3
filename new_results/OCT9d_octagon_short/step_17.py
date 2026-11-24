import os, sys, pickle, time
base = "<project>/research/08_theory/"
O9B = base+"OCT9b_octagon_debug"
sys.path.insert(0, O9B); sys.path.insert(0, base+"OCT9c_octagon_close")
import sympy as sp
from fractions import Fraction as Fr
from syminfra import S_Li3, S_mLi2logp, S_logprod, SymbolSpace, symbol_of_expr
from oct9c_common import load_all, get_engine, rhs_terms_235, rhs_terms_2357

D = load_all()
# Symbol spaces
SS235 = SymbolSpace((2,3,5))     # dim 27
SS2357 = SymbolSpace((2,3,5,7))  # dim 64
print("dims:", SS235.dim3, SS2357.dim3)

t0=time.time()
TW3 = {}   # target symbol vectors
for e in range(4):
    eng = get_engine(e)
    if e == 0:
        W3 = D['Wp0']['W3_e0']; SS = SS2357
    else:
        W3 = D['s123'][e]['W3']; SS = SS235
    S = symbol_of_expr(W3, eng.ATOMS)
    v = SS.coord(S)
    TW3[e] = v
    nz = sum(1 for x in v if x != 0)
    print(f"edge {e}: symbol nonzero coords = {nz} / {SS.dim3}   ({time.time()-t0:.1f}s)")
pickle.dump({e: [sp.Rational(x) for x in TW3[e]] for e in TW3}, open("TW3.pkl","wb"))
print("saved TW3.pkl")
