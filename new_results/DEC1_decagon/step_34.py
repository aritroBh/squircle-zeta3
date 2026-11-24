import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 80
s5 = mp.sqrt(5); phi = (1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
s36 = mp.sqrt(10-2*s5)/4; R=1/(2*phi**2*s36); h=1/(2*phi); halfL=R/(2*phi); L=2*halfL
E3n = pickle.load(open('E3n_num100.pkl','rb'))
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
invn=[mp.mpf(1),mp.mpf(-1),mp.mpf(1),mp.mpf(-1)]
m2lg2=[-2*x for x in lg2]
fac2=npexp(npmul(m2lg2,invn))
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
    f.write(f"step: E3(0.3) check vs direct d-series: diff = {mp.nstr(abs(E3_engine-E3_direct),4)} (80 dps)\n")