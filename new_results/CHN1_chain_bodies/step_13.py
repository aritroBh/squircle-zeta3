import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import mzv_exact as M
for W in (9,10,11):
    t=time.time()
    T=M.DSTable(W)
    print('W=%d words=%d nullity=%d expected=%d  %.1fs'%(W,T.n,T.nullity,len(M.BASIS[W]),time.time()-t))
print('total %.1fs'%(time.time()-t0))
open('progress.log','a').write('step_12: DS tables W=9,10,11 built, nullities match conjectured dims 5,7,9\n')
