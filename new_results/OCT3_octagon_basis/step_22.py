import json, time
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr, nint
mp.dps = 280
from flint import fmpz_mat
import sympy as sp
t0=time.time()
r2=sqrt(2); b=r2-1; u=sqrt(2+r2)
l2=log(2); lu=log(1+r2); z3=zeta(3); pi2=pi**2
Pk=[mpf(s) for s in json.load(open("pder.json"))["vals"]]
zs=[-1,-r2/2,-b,(1-r2)/2,-b**2,(4-3*r2)/8,b**4,(2-r2)/4,b**2,1-r2/2,b,mpf('0.5'),2-r2,r2/2,2*b,(2+r2)/4,12*r2-16]
n2s=["-1","-r2/2","-b","(1-r2)/2","-b^2","(4-3r2)/8","b^4","(2-r2)/4","b^2","1-r2/2","b","1/2","2-r2","r2/2","2b","(2+r2)/4","12r2-16"]
Li2=[polylog(2,z) for z in zs]
w2=Li2+[pi2,l2**2,l2*lu,lu**2]
w2n=["Li2(%s)"%s for s in n2s]+["pi2","l2^2","l2*lu","lu^2"]
X2=[];N2=[]
for nm,c in zip(w2n,w2): X2+=[c,r2*c]; N2+=[nm,"r2*"+nm]
n=len(X2);B=10**255
rows=[]
for i in range(n):
    r=[0]*(n+1);r[i]=1;r[n]=int(nint(B*X2[i]));rows.append(r)
Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
rel=[]
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,X2) if ai)
    if abs(s)<mpf('1e-270'): rel.append(a)
A=sp.Matrix(rel); R,piv=A.rref()
indep=[j for j in range(n) if j not in piv]
print("independent weight-2 constants:")
for j in indep: print("  ",N2[j])
# expressions for a couple of eliminated Li2
for i,p in enumerate(piv):
    if "Li2" in N2[p] and len([1 for j in indep if R[i,j]!=0])>0:
        terms=[(R[i,j],N2[j]) for j in indep if R[i,j]!=0]
        if N2[p] in ("Li2(b)","Li2(1/2)","Li2(-1)"):
            print(N2[p],"="," ".join("%s*%s"%(-c,nm) for c,nm in terms))

# hunt P'' : basis = (indep weight-2 base) + {l2,lu,1} x 4 mults
w2base = [X2[j] for j in indep if not N2[j].startswith("r2")]
w2basen = [N2[j] for j in indep if not N2[j].startswith("r2")]
low = [l2,lu,mpf(1)]; lown=["l2","lu","1"]
base = w2base+low; basen = w2basen+lown
mults=[(mpf(1),"1"),(r2,"r2"),(u,"u"),(u*r2,"r2*u")]
XB=[];NB=[]
for nm,c in zip(basen,base):
    for mv,mn in mults:
        XB.append(c*mv); NB.append(nm if mn=="1" else mn+"*"+nm)
print("P'' basis size:", len(XB))
def hunt(target,tname,vals,vnames,B=10**260,top=6):
    v=[target]+vals;n=len(v)
    rows=[]
    for i in range(n):
        r=[0]*(n+1);r[i]=1;r[n]=int(nint(B*v[i]));rows.append(r)
    Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
    cs=[]
    for row in Rl:
        a=row[:n]
        s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
        cs.append((abs(s),max(abs(ai) for ai in a),abs(a[0]),a))
    cs.sort(key=lambda t:t[0])
    print("hunt %s (dim %d):"%(tname,n))
    for c in cs[:top]: print("   |s|=%s maxc=%-10d |a0|=%d"%(nstr(c[0],3),c[1],c[2]))
    return cs
cs=hunt(Pk[2],"P'' vs weight2-basis",XB,NB)
print("[%.1fs]"%(time.time()-t0))
