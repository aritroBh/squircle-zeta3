import os, pickle; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 280
d=pickle.load(open('basis280.pkl','rb')); indep=pickle.load(open('indep280.pkl','rb'))
allv=[mp.mpf(x) for x in d['allv']]; alln=d['alln']
r2=mp.sqrt(2)
bvals=[allv[i] for i in indep]; bn=[alln[i] for i in indep]
D1=mp.mpf(d['D1']); D2=mp.mpf(d['D2']); P3u=mp.mpf(d['P3u']); P3=mp.mpf(d['P3'])
from flint import fmpz_mat
def hunt(target, consts, cn, dps, se, withsqrt2=True):
    old=mp.mp.dps; mp.mp.dps=dps; B=mp.power(10,se)
    vs=[target]; nm=['T']
    for c,nme in zip(consts,cn):
        vs.append(c); nm.append(nme)
        if withsqrt2: vs.append(c*r2); nm.append(nme+'r2')
    n=len(vs)
    rows=[[0]*(n+1) for _ in range(n)]
    for i in range(n): rows[i][i]=1; rows[i][-1]=int(mp.nint(B*vs[i]))
    R=fmpz_mat(rows).lll(); out=[]
    for i in range(R.nrows()):
        v=[int(R[i,j]) for j in range(R.ncols())]
        if v[0]!=0: out.append((abs(v[-1]),max(abs(cc) for cc in v[:-1]),v))
    mp.mp.dps=old; out.sort(); return out, nm
# Hunt 1: P3u over {1,r2} x indep singles (no D1,D2)
res,nm = hunt(P3u, bvals, bn, 280, 250)
print("P3u over {1,r2}x20indep (no D1,D2): floor maxc=%.2e residB=%.1e"%(res[0][1],res[0][0]))
# Hunt 2: add D1,D2
res2,nm2 = hunt(P3u, bvals+[D1,D2], bn+['D1','D2'], 280, 250)
r0=res2[0]
print("P3u over {1,r2}x(20indep+D1+D2): floor maxc=%.2e residB=%.1e"%(r0[1],r0[0]))
if r0[0]<1000:
    combo=[(nm2[i],r0[2][i]) for i in range(1,len(r0[2])-1) if r0[2][i]!=0]
    print("  RELATION:", combo)
flog("step: hunted P3u at 280dps with/without D1,D2")
