import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import mzv_exact as M
for W in range(2,9):
    t=time.time()
    T=M.DSTable(W)
    print('W=%d words=%d relations_nullity=%d expected=%d  %.1fs'%(W,T.n,T.nullity,len(M.BASIS[W]),time.time()-t))
print('total %.1fs'%(time.time()-t0))
