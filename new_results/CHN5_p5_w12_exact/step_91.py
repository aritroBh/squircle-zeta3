import os, sys, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
exec(open('run_gamma_trace2.py').read().split('schedule=')[0])  # load solver defs (QA,QC,v10_mpf,gamma_at_mpf)
t0=time.time()
# probe large-a1 end
for a1v in [3e6, 1e7, 3e7]:
    sols=gamma_at_mpf(a1v, [], n=60, lo=-4, hi=7)
    print(f"a1={a1v:.3g}: {len(sols)} pts", [('+' if s['QA']>0 else '-', '+' if s['QC']>0 else '-') for s in sols])
    for s in sols:
        print(f"   a2={mp.nstr(s['a2'],8)} a4={mp.nstr(s['a4'],8)} a5={mp.nstr(s['a5'],8)} QA/S^8={mp.nstr(s['QA']/(s['a1']+s['a2']+1+s['a4']+s['a5'])**8,4)}")
print(f"[{time.time()-t0:.0f}s]")
