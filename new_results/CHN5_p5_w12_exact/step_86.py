import os, pickle, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
mp.mp.dps=60
res=pickle.load(open('gamma_trace2_partial.pkl','rb'))
for k in res:
    if abs(k-1.0)/1.0 < 1e-9:
        for s in res[k]:
            print(f"a1={mp.nstr(s['a1'],8)} a2={mp.nstr(s['a2'],12)} a4={mp.nstr(s['a4'],12)} a5={mp.nstr(s['a5'],12)}")
            S=s['a1']+s['a2']+1+s['a4']+s['a5']
            print(f"   QA={mp.nstr(s['QA'],8)}  QC={mp.nstr(s['QC'],8)}  QA/S^8={mp.nstr(s['QA']/S**8,4)} QC/S^8={mp.nstr(s['QC']/S**8,4)}")
