import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import mzv_exact as M, engine as E
def table(W):
    if W not in E._TABLES: E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
for a in [(1,1,1,1),(3,1,2,5)]:
    r=E.reduced_coeffs(a,11)
    print(a,'W10 z2z35=%s z37=%s | W11 z3z35=%s z335=%s'%(r[10][5],r[10][6],r[11][7],r[11][8]))
print('t=%.1fs'%(time.time()-t0))
