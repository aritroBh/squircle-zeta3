"""function_symbol_check.py -- EXACT.
The function-level symbol of a polylogarithmic function of c is the constant-retaining symbol with every
constant letter (here: 2) set to 0 (constants such as log 2 are then complex scalars, not weight-1 letters).
This is the version needed for an argument about FUNCTIONS of c: equal functions have equal function symbols,
classical Li_4(f(c)), products and constants x lower-weight functions contribute 0 to delta_{2,2} of the
weight-4 part.  We recompute delta_{2,2} of the three prefactor components of W_4(c) in this version, and the
rank of the obstructions of A_{-1}(c), A_0(c), A_{-1}(s), A_0(s) after pullback to tau (constant letter dropped).
"""
import sys, itertools
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, '.')
import hpl
from b2b2_test import load, lin, delta22
import obstruction_structure as O

def drop_const(S):
    return {k: v for k, v in S.items() if 0 not in k}   # index 0 = letter '2'

for e in (-1, 0, 1):
    D = delta22(drop_const(O.comp[e]))
    print(f'function-symbol delta22(A_{e}) : {len(D)} nonzero components: {D}')

# tau pullback with the constant letter (index 0 = 2) removed
def drop5(S):
    return {k: v for k, v in S.items() if 0 not in k}
obs = {}
for e in (-1, 0):
    for nm, xf in (('c', O.cfun), ('s', O.sfun)):
        obs[(e, nm)] = O.delta22_5(drop5(O.pullback(drop_const(O.comp[e]), xf)))
keys = [(-1, 'c'), (0, 'c'), (-1, 's'), (0, 's')]
rows = sorted(set().union(*[set(obs[k]) for k in keys]))
Mx = sp.Matrix([[obs[k].get(r, 0) for k in keys] for r in rows])
print('function-symbol obstruction rank of A_{-1}(c), A_0(c), A_{-1}(s), A_0(s) in the tau alphabet:', Mx.rank())
for e in (1,):
    for nm, xf in (('c', O.cfun), ('s', O.sfun)):
        print(f'function-symbol delta22(A_1({nm}(tau))) components:', len(O.delta22_5(drop5(O.pullback(drop_const(O.comp[e]), xf)))))
