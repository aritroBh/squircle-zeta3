import os, sys, pickle, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
mp.mp.dps=150
res=pickle.load(open('gamma_trace2.pkl','rb'))
for k in [0.01, 1.7782794100389228, 316.22776601683796]:
    for s in res[k]:
        if s['a5']>mp.mpf('1e-25'):
            print(f"a1={k}: a2={mp.nstr(s['a2'],12)} a4={mp.nstr(s['a4'],12)} a5={mp.nstr(s['a5'],12)}")
