import json, time
from mpmath import mp, mpf, sqrt, log, pi, zeta, polylog, nstr, nint
mp.dps = 280
from flint import fmpz_mat
t0=time.time()
r2=sqrt(2); b=r2-1; u=sqrt(2+r2)
l2=log(2); lu=log(1+r2); z3=zeta(3); pi2=pi**2
Pk=[mpf(s) for s in json.load(open("pder.json"))["vals"]]

def hunt(target, tname, vals, vnames, B=10**260, top=6):
    v=[target]+vals; n=len(v)
    rows=[]
    for i in range(n):
        r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*v[i])); rows.append(r)
    Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
    cs=[]
    for row in Rl:
        a=row[:n]
        s=sum(mpf(ai)*xi for ai,xi in zip(a,v) if ai)
        cs.append((abs(s), max(abs(ai) for ai in a), abs(a[0]), a, s))
    cs.sort(key=lambda t:t[0])
    print("hunt %s (dim %d):"%(tname,n))
    for c in cs[:top]:
        print("   |s|=%s maxc=%-10d |a0|=%d"%(nstr(c[0],3),c[1],c[2]))
    return cs

# A. lean weight-3 basis x {1,r2,u,r2u}
Li3b = polylog(3,b); Li3h = polylog(3,r2/2)
leanbase = [Li3b, Li3h, z3, pi2*l2, pi2*lu, l2**3, l2**2*lu, l2*lu**2, lu**3,
            pi2, l2**2, l2*lu, lu**2, l2, lu, mpf(1)]
leannames = ["Li3(b)","Li3(r2/2)","z3","pi2*l2","pi2*lu","l2^3","l2^2*lu","l2*lu^2","lu^3",
             "pi2","l2^2","l2*lu","lu^2","l2","lu","1"]
mults=[(mpf(1),"1"),(r2,"r2"),(u,"u"),(u*r2,"r2*u")]
XL=[]; NL=[]
for nm,c in zip(leannames,leanbase):
    for mv,mn in mults:
        XL.append(c*mv); NL.append(nm if mn=="1" else mn+"*"+nm)
cs = hunt(Pk[3], "P''' vs lean64", XL, NL)
print("  [%.1fs]"%(time.time()-t0))

# B. weight-2 layer
Li2 = [polylog(2,z) for z in [-1,-r2/2,-b,(1-r2)/2,-b**2,(4-3*r2)/8,b**4,(2-r2)/4,
                              b**2,1-r2/2,b,mpf('0.5'),2-r2,r2/2,2*b,(2+r2)/4,12*r2-16]]
n2=["-1","-r2/2","-b","(1-r2)/2","-b^2","(4-3r2)/8","b^4","(2-r2)/4","b^2","1-r2/2","b","1/2","2-r2","r2/2","2b","(2+r2)/4","12r2-16"]
w2 = Li2 + [pi2, l2**2, l2*lu, lu**2]
w2n = ["Li2(%s)"%s for s in n2]+["pi2","l2^2","l2*lu","lu^2"]
X2=[]; N2=[]
for nm,c in zip(w2n,w2):
    X2 += [c, r2*c]; N2 += [nm, "r2*"+nm]
# internal relations on weight-2 (42 constants)
rows=[]; n=len(X2); B=10**255
for i in range(n):
    r=[0]*(n+1); r[i]=1; r[n]=int(nint(B*X2[i])); rows.append(r)
Rl=[[int(x) for x in row] for row in fmpz_mat(rows).lll().tolist()]
rel=[]
for row in Rl:
    a=row[:n]
    s=sum(mpf(ai)*xi for ai,xi in zip(a,X2) if ai)
    if abs(s)<mpf('1e-270'): rel.append((a,s))
print("\nweight-2 internal relations found:", len(rel), " -> independent:", n-len(rel))
for a,s in rel[:30]:
    print("  RES",nstr(s,1),"|"," ".join("%+d*%s"%(ai,N2[i]) for i,ai in enumerate(a) if ai))
print("  [%.1fs]"%(time.time()-t0))
