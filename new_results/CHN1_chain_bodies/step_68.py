import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import mzv_exact as M, engine as E
from fractions import Fraction as F
# patch engine.table to use flint for W>=11
import types
_orig=E.table
def table(W):
    if W not in E._TABLES:
        E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
t=time.time()
r=E.reduced_coeffs((2,2,3,3),12)[12]
print('one W=12 engine call: %.0fs'%(time.time()-t))
labels=['z12','z3z9','z2z3z7','z2z5^2','z2^2z3z5','z2^3z3^2','z3^4','z2^2z35','z2z37','z39','z57','w*']
print('Pi1 point (2,2,3,3), c12 non-product coords:')
for i in (9,10,11):
    print(' ', labels[i], '=', r[i])
print('  irreducible-products:', labels[7],'=',r[7], labels[8],'=',r[8])
print('t=%.1fs'%(time.time()-t0))
