import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=340
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
d=pickle.load(open(WORK+'/cands29.pickle','rb'))
fw=d['freewords']; allv=d['cands']
out=pickle.load(open(WORK+'/vals110.pickle','rb'))
V={k:[mp.mpc(mp.mpf(out['vals'][k][w][0]),mp.mpf(out['vals'][k][w][1])) for w in fw] for k in (1,3,5,7)}
def res(row,k):
    s=mp.mpc(0)
    for i,c in enumerate(row):
        if c: s+=c*V[k][i]
    return abs(s)
true29=[]; acc=[]
for row in allv:
    r=res(row,1)
    lr=float(mp.log10(r)) if r>0 else -999
    (true29 if lr<-320 else acc).append(row)
print('true (resid<1e-320):',len(true29),' accidents:',len(acc))
print('accident resid range:', f"{min(float(mp.log10(res(r,1))) for r in acc):.1f} .. {max(float(mp.log10(res(r,1))) for r in acc):.1f}")
print('true resid range:', f"{max(float(mp.log10(res(r,1))) for r in true29):.1f}")
# classify by k=3,5
print('idx supp height   k=1     k=3     k=5     k=7')
nform=0
for i,row in enumerate(true29):
    rs=[res(row,k) for k in (1,3,5,7)]
    lrs=[float(mp.log10(r)) if r>0 else -999 for r in rs]
    form = lrs[1]<-300 and lrs[2]<-300
    if form: nform+=1
    print(f'{i:3d} {sum(1 for x in row if x):4d} {max(abs(x) for x in row):7.1e} '+' '.join(f'{x:7.1f}' for x in lrs)+('   FORM' if form else ''))
print('form-level count:',nform)
pickle.dump({'freewords':fw,'cands':true29},open(WORK+'/cands29.pickle','wb'))
print('saved 29 candidates')
