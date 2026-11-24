import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import mzv_exact as M, engine as E
def table(W):
    if W not in E._TABLES:
        E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
labels=['z12','z3z9','z2z3z7','z2z5^2','z2^2z3z5','z2^3z3^2','z3^4','z2^2z35','z2z37','z39','z57','w*']
for a in [(2,3,2,3),(1,1,1,1),(1,2,3,4),(3,1,1,1),(2,1,1,3)]:
    t=time.time()
    r=E.reduced_coeffs(a,12)[12]
    loc = 'Pi1' if a[0]==a[1] and a[2]==a[3] else ('Pi2' if a[0]==a[2] and a[1]==a[3] else ('diag' if a[0]==a[1]==a[2]==a[3] else 'generic'))
    print(a, loc, 'z39=%s z57=%s w*=%s z2^2z35=%s z2z37=%s (%.0fs)'%(r[9],r[10],r[11],r[7],r[8],time.time()-t))
print('t=%.1fs'%(time.time()-t0))
