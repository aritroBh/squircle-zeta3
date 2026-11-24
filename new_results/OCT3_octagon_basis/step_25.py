import json, time
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr, nint
mp.dps = 280
from flint import fmpz_mat
import sympy as sp
t0=time.time()
r2=sqrt(2); b=r2-1; u=sqrt(2+r2)
l2=log(2); lu=log(1+r2); z3=zeta(3); pi2=pi**2
Pk=[mpf(s) for s in json.load(open("pder.json"))["vals"]]
Li2g=polylog(2,12*r2-16)
Li3_2b=polylog(3,2*b); Li3_q=polylog(3,(2+r2)/4); Li3_t=polylog(3,12*r2-16)

w3=[Li3_2b,Li3_q,Li3_t,z3,pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]
w3n=["Li3(2b)","Li3((2+r2)/4)","Li3(12r2-16)","z3","pi2*l2","pi2*lu","l2^3","l2^2*lu","l2*lu^2","lu^3"]
w2=[pi2,l2**2,l2*lu,lu**2,Li2g]
w2n=["pi2","l2^2","l2*lu","lu^2","Li2(12r2-16)"]
w10=[l2,lu,mpf(1)]; w10n=["l2","lu","1"]
base=w3+w2+w10; basen=w3n+w2n+w10n
XB=[];NB=[]
for nm,c in zip(basen,base):
    XB += [c, r2*c]; NB += [nm, "r2*"+nm]
print("basis size:", len(XB))

def hunt(target,tname,vals,vnames,B=10**260,top=5):
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

T = Pk[3]/u
cs = hunt(T, "P'''/u vs 36", XB, NB)
if cs[0][0] < mpf('1e-255'):
    a=cs[0][3]
    print("\nFOUND relation, residual", nstr(cs[0][0],3))
    vn=["P'''/u"]+NB
    a0=sp.Integer(a[0])
    print("a0 =", a0)
    for i,ai in enumerate(a):
        if ai: print("   %+8d * %s"%(ai,vn[i]))
print("[%.1fs]"%(time.time()-t0))
