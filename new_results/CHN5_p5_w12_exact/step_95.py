import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])
mp.mp.dps=70
def q11_signs_present(a1v, lo, hi, n=80):
    sgns=set(); cnt=0
    for a2v in np.logspace(lo,hi,n):
        for (a4v,a5v) in v10_mpf(a1v,a2v):
            q=q11_mpf([a1v,a2v,1,a4v,a5v]); sgns.add('+' if q>0 else '-'); cnt+=1
    return sgns, cnt
t0=time.time()
for a1v in [1e8, 1e9, 1e10]:
    lo=math.log10(a1v)-3 if (math:=__import__('math')) else 0
    sgns,cnt=q11_signs_present(a1v, lo, lo+3)
    print(f"a1={a1v:.0e}: V10 sheet pts={cnt}, Q11 signs present={sgns} [{time.time()-t0:.0f}s]")
