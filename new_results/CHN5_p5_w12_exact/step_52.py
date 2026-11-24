import os
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
analyze = r'''
import pickle, mpmath as mp
mp.mp.dps=60
res=pickle.load(open('gamma_trace.pkl','rb'))
a1s=sorted(res.keys())
# collect all solution points
allpts=[]
for a1v in a1s:
    for s in res[a1v]:
        allpts.append(s)
print(f"total Gamma points traced: {len(allpts)} over {len(a1s)} a1 values [{mp.nstr(min(a1s),3)} .. {mp.nstr(max(a1s),3)}]")
# link into branches: sort by a1, nearest-neighbor in log(a2),log(a4),log(a5)
import math
def key(s):
    return (math.log(float(s['a2'])),math.log(float(s['a4'])),math.log(float(s['a5'])))
branches=[]
for a1v in a1s:
    for s in res[a1v]:
        k=key(s)
        best=None; bestd=1e9
        for br in branches:
            lk=br['last']
            d=sum((k[i]-lk[0][i])**2 for i in range(3))
            if d<bestd: bestd=d; best=br
        if best is not None and bestd< (0.35*3)**2 and best['lasta1']<a1v*1.3:
            best['pts'].append(s); best['last']=(k,); best['lasta1']=a1v
        else:
            branches.append({'pts':[s],'last':(k,),'lasta1':a1v})
big=[b for b in branches if len(b['pts'])>=3]
print(f"branches (>=3 pts): {len(big)}")
minnormQA=mp.mpf('1e100'); minarg=None
for i,b in enumerate(big):
    pts=b['pts']
    sgnQA=set('+' if p['QA']>0 else '-' for p in pts)
    sgnQC=set('+' if p['QC']>0 else '-' for p in pts)
    a1lo=float(pts[0]['a1']); a1hi=float(pts[-1]['a1'])
    # normalized |QA| by (sum a)^8
    for p in pts:
        S=float(p['a1']+p['a2']+1+p['a4']+p['a5'])
        n=abs(float(p['QA']))/S**8
        if n<minnormQA: minnormQA=n; minarg=(i,p)
    print(f" branch {i}: {len(pts)} pts, a1 in [{a1lo:.3g},{a1hi:.3g}], signQA={sorted(sgnQA)}, signQC={sorted(sgnQC)}")
    # check sign changes
    for j in range(1,len(pts)):
        if pts[j]['QA']*pts[j-1]['QA']<0:
            print(f"    *** QA SIGN CHANGE between a1={mp.nstr(pts[j-1]['a1'],6)} and {mp.nstr(pts[j]['a1'],6)}")
        if pts[j]['QC']*pts[j-1]['QC']<0:
            print(f"    *** QC SIGN CHANGE between a1={mp.nstr(pts[j-1]['a1'],6)} and {mp.nstr(pts[j]['a1'],6)}")
print(f"min normalized |QA|/(S^8) = {mp.nstr(mpnorm if (mpnorm:=minnormQA) else 0,6)}")
if minarg:
    i,p=minarg
    print(f"  at branch {i}: a1={mp.nstr(p['a1'],10)} a2={mp.nstr(p['a2'],10)} a4={mp.nstr(p['a4'],10)} a5={mp.nstr(p['a5'],10)}")
    print(f"  QA={mp.nstr(p['QA'],10)} QC={mp.nstr(p['QC'],10)}")
'''
open(W+'/run_gamma_analyze.py','w').write(analyze)
print("run_gamma_analyze.py written")
