import os, pickle, mpmath as mp, math
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
mp.mp.dps=60
res=pickle.load(open('gamma_trace2.pkl','rb'))
a1s=sorted(res.keys())
# filter stratum artifacts: keep points with a5 > 1e-25 (genuine interior); also require not (a2==a4 exact & a5<1e-20)
tot=0; kept=0; dropped=0
minnormQA=mp.mpf('1e100'); minarg=None
minnormQC=mp.mpf('1e100'); minargC=None
signsQA=set(); signsQC=set()
for a1v in a1s:
    for s in res[a1v]:
        tot+=1
        if s['a5']<mp.mpf('1e-25'):
            dropped+=1; continue
        kept+=1
        S=s['a1']+s['a2']+1+s['a4']+s['a5']
        nA=s['QA']/S**8; nC=s['QC']/S**8
        signsQA.add('+' if s['QA']>0 else '-'); signsQC.add('+' if s['QC']>0 else '-')
        if abs(nA)<minnormQA: minnormQA=abs(nA); minarg=(a1v,s,nA)
        if abs(nC)<minnormQC: minnormQC=abs(nC); minargC=(a1v,s,nC)
print(f"total pts {tot}, kept {kept}, dropped stratum {dropped}")
print(f"signQA (interior, filtered): {sorted(signsQA)}, signQC: {sorted(signsQC)}")
print(f"min |QA|/S^8 = {mp.nstr(minnormQA,6)} at a1={mp.nstr(minarg[0],8)} a2={mp.nstr(minarg[1]['a2'],8)} a5={mp.nstr(minarg[1]['a5'],6)} (QA={mp.nstr(minarg[1]['QA'],5)})")
print(f"min |QC|/S^8 = {mp.nstr(minnormQC,6)} at a1={mp.nstr(minargC[0],8)} (QC={mp.nstr(minargC[1]['QC'],5)})")
# branch structure: count interior pts per a1
print("\ninterior branch counts (filtered):")
import numpy as np
counts=[(a1v,len([s for s in res[a1v] if s['a5']>mp.mpf('1e-25')])) for a1v in a1s]
print(" ".join(f"{c}" for _,c in counts))
