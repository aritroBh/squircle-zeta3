import os, mpmath as mp, pickle, math
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
s5 = mp.sqrt(5); phi = (1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
h2oL2=(5+2*s5)/4
def eadd(P,Q):
    R={}
    for k in set(P)|set(Q):
        aa=P.get(k,[mp.mpf(0)]*4); b=Q.get(k,[mp.mpf(0)]*4)
        R[k]=[aa[i]+b[i] for i in range(4)]
    return R
def ekmul(P,Q):
    R={}
    for k1,v1 in P.items():
        for k2,v2 in Q.items():
            k=tuple(k1[i]+k2[i] for i in range(5))
            if sum(k)>3: continue
            c=[mp.mpf(0)]*4
            for i in range(4):
                for j in range(4-i): c[i+j]+=v1[i]*v2[j]
            if k in R: R[k]=[R[k][i]+c[i] for i in range(4)]
            else: R[k]=c
    return R
def escale(P,s): return {k:[v[i]*s for i in range(4)] for k,v in P.items()}
def ediv_scalar(P,s): return {k:[v[i]/s for i in range(4)] for k,v in P.items()}
ONE={(0,0,0,0,0):[mp.mpf(1),0,0,0]}
def pmul(P,Q):
    R=[{} for _ in range(4)]
    for i in range(4):
        for j in range(4-i): R[i+j]=eadd(R[i+j],ekmul(P[i],Q[j]))
    return R
def pdiv(P,Q):
    R=[{} for _ in range(4)]
    q0=Q[0][(0,0,0,0,0)][0]
    for i in range(4):
        acc=P[i]
        for j in range(i): acc=eadd(acc,escale(ekmul(R[j],Q[i-j]),mp.mpf(-1)))
        R[i]=ediv_scalar(acc,q0)
    return R
def plog(P):
    Pp=[escale(P[i],mp.mpf(i)) for i in range(1,4)]+[{}]
    Q=pdiv(Pp,P)
    L=[{}]
    for i in range(3): L.append(escale(Q[i],mp.mpf(1)/(i+1)))
    return L
def pexp(P):
    E=[ONE]+[{}]*3
    for k in range(1,4):
        acc={}
        for j in range(1,k+1): acc=eadd(acc,escale(ekmul(P[j],E[k-j]),mp.mpf(j)))
        E[k]=escale(acc,mp.mpf(1)/k)
    return E
apol = [[mp.mpf(0),A0],[A1,-A1],[B1,B1*kap1],[B4,B4*kap4],[B5,mp.mpf(0)]]
appol = [[A0,mp.mpf(0)],[-A1,mp.mpf(0)],[B1*kap1,mp.mpf(0)],[B4*kap4,mp.mpf(0)],[mp.mpf(0),mp.mpf(0)]]
F=[{} for _ in range(4)]; S=[{} for _ in range(4)]
for k in range(4):
    fc=mp.mpf(1)/math.factorial(k)
    for i in range(5):
        key=tuple(k if j==i else 0 for j in range(5))
        F[k]=eadd(F[k],{key:[apol[i][0]*fc, apol[i][1]*fc,0,0]})
        S[k]=eadd(S[k],{key:[appol[i][0]*fc, appol[i][1]*fc,0,0]})
Phi=pdiv(S,F); lg=plog(F)
inv=[escale(ONE,mp.mpf(s)) for s in [1,-1,1,-1]]
m2lg=[escale(x,mp.mpf(-2)) for x in lg]
facE=pexp(pmul(m2lg,inv))
Phisq=pmul(Phi,Phi)
Q=[{} for _ in range(4)]
Q[0]=ONE
o2mt = {(0,0,0,0,0):[mp.mpf(1),mp.mpf(-2),0,0]}   # (1-2t)
yy2poly=[h2oL2+mp.mpf('0.25'), mp.mpf(-1), mp.mpf(1), mp.mpf(0)]
for k in range(1,4): Q[k]=ekmul(Phi[k], o2mt)      # -2*yypn(t)*Phi[k], yypn=-(1-2t)/2
for k in range(4):
    Q[k]=eadd(Q[k], ekmul(Phisq[k], {(0,0,0,0,0):yy2poly}))
G=pmul(facE,Q)
sq=[{} for _ in range(4)]; sq[0]=ONE
for k in range(1,4):
    acc=G[k]
    for j in range(1,k): acc=eadd(acc,escale(ekmul(sq[j],sq[k-j]),mp.mpf(-1)))
    sq[k]=escale(acc,mp.mpf('0.5'))
E3n=sq[3]
pickle.dump(E3n, open('E3n_num100.pkl','wb'))
# ---- check at t=0.3 ----
mp.mp.dps=80
s36 = mp.sqrt(10-2*s5)/4; R=1/(2*phi**2*s36); h=1/(2*phi); halfL=R/(2*phi); L=2*halfL
t0 = mp.mpf('0.3')
avals = [A0*t0, A1*(1-t0), B1*(1+kap1*t0), B4*(1+kap4*t0), B5]
logs = [mp.log(v) for v in avals]
val = mp.mpf(0)
for key,poly in E3n.items():
    mval = poly[0]+poly[1]*t0+poly[2]*t0**2+poly[3]*t0**3
    lv = mp.mpf(1)
    for i,e in enumerate(key): lv *= logs[i]**e
    val += mval*lv
E3_engine = L*val
y = mp.matrix([-h, halfL*(1-2*t0)]); yp = mp.matrix([mp.mpf(0), -2*halfL])
a_t = avals; ap_t = [A0, -A1, B1*kap1, B4*kap4, mp.mpf(0)]
n=3
def npmul(P,Q):
    R=[mp.mpf(0)]*(n+1)
    for i in range(n+1):
        for j in range(n+1-i): R[i+j]+=P[i]*Q[j]
    return R
def npdiv(P,Q):
    R=[mp.mpf(0)]*(n+1)
    for i in range(n+1): R[i]=(P[i]-sum(R[j]*Q[i-j] for j in range(i)))/Q[0]
    return R
def nplog(P):
    Pp=[P[i]*i for i in range(1,n+1)]+[mp.mpf(0)]
    Q=npdiv(Pp,P); Ls=[mp.mpf(0)]
    for i in range(n): Ls.append(Q[i]/(i+1))
    return Ls
def npexp(P):
    E=[mp.mpf(1)]+[mp.mpf(0)]*n
    for k in range(1,n+1): E[k]=sum(P[j]*j*E[k-j] for j in range(1,k+1))/k
    return E
Fs=[mp.mpf(0)]*(n+1); Ss=[mp.mpf(0)]*(n+1)
for aj,lj in zip(a_t,logs):
    pk=mp.mpf(1)
    for k in range(n+1):
        if k: pk=pk*lj/k
        Fs[k]+=aj*pk
for aj,lj in zip(ap_t,logs):
    pk=mp.mpf(1)
    for k in range(n+1):
        if k: pk=pk*lj/k
        Ss[k]+=aj*pk
Phi2=npdiv(Ss,Fs); lg2=nplog(Fs)
m2lg2=[-2*x for x in lg2]
fac2=npexp(npmul(m2lg2,inv))
yp2=yp[0]**2+yp[1]**2; yyp=y[0]*yp[0]+y[1]*yp[1]; yy2=y[0]**2+y[1]**2
Phisq2=npmul(Phi2,Phi2)
Q2=[mp.mpf(0)]*(n+1); Q2[0]=yp2
for k in range(1,n+1): Q2[k]=-2*yyp*Phi2[k]
for k in range(n+1): Q2[k]+=yy2*Phisq2[k]
G2=npmul(fac2,Q2)
sq2=[mp.mpf(0)]*(n+1); sq2[0]=mp.sqrt(G2[0])
for k in range(1,n+1): sq2[k]=(G2[k]-sum(sq2[j]*sq2[k-j] for j in range(1,k)))/(2*sq2[0])
E3_direct = sq2[3]
print("E3(0.3) engine:", mp.nstr(E3_engine,40))
print("E3(0.3) direct:", mp.nstr(E3_direct,40))
print("diff:", mp.nstr(abs(E3_engine-E3_direct),6))
with open('progress.log','a') as f:
    f.write(f"step: E3(0.3) check vs direct d-series (yyp(t) fixed): diff = {mp.nstr(abs(E3_engine-E3_direct),4)} (80 dps)\n")