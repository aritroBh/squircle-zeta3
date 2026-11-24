import os, pickle, mpmath as mp, math
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
mp.mp.dps=60
res=pickle.load(open('gamma_trace2.pkl','rb'))
a1s=sorted(res.keys())
# collect interior points with normalized QA,QC; link into branches by nearest neighbor
pts=[]
for a1v in a1s:
    for s in res[a1v]:
        if s['a5']>mp.mpf('1e-25'):
            S=s['a1']+s['a2']+1+s['a4']+s['a5']
            pts.append({'a1':a1v,'a2':s['a2'],'a4':s['a4'],'a5':s['a5'],'nQA':s['QA']/S**8,'nQC':s['QC']/S**8})
def key(p): return (math.log(float(p['a2'])),math.log(float(p['a4'])),math.log(float(p['a5'])))
# simple branch linking sorted by a1
branches=[]
for p in pts:
    k=key(p); best=None; bd=9e9
    for br in branches:
        lk=br[-1]
        d=sum((k[i]-lk['k'][i])**2 for i in range(3))
        if d<bd: bd=d; best=br
    if best is not None and bd<4.0:
        best.append({'k':k,'p':p})
    else:
        branches.append([{'k':k,'p':p}])
branches=[b for b in branches if len(b)>=2]
print(f"{len(pts)} interior pts, {len(branches)} branches (>=2 pts)")
for bi,b in enumerate(branches):
    ps=[e['p'] for e in b]
    a1lo,a1hi=ps[0]['a1'],ps[-1]['a1']
    nQA=[p['nQA'] for p in ps]
    # interior local minima of |nQA|
    dips=[]
    for i in range(1,len(ps)-1):
        if abs(nQA[i])<abs(nQA[i-1]) and abs(nQA[i])<abs(nQA[i+1]):
            dips.append((ps[i]['a1'],nQA[i]))
    sgn=set('+' if q>0 else '-' for q in nQA)
    print(f"br{bi}: {len(ps)} pts a1[{mp.nstr(a1lo,3)},{mp.nstr(a1hi,3)}] sign nQA={sgn} min|nQA|={mp.nstr(min(abs(q) for q in nQA),3)}" + (f" DIPS={[(mp.nstr(a,3),mp.nstr(v,3)) for a,v in dips]}" if dips else ""))
