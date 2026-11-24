import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 80
s5 = mp.sqrt(5); phi = (1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
s36 = mp.sqrt(10-2*s5)/4; R=1/(2*phi**2*s36); h=1/(2*phi)
halfL=R/(2*phi); L=2*halfL
E3n = pickle.load(open('E3n_num100.pkl','rb'))
t0 = mp.mpf('0.3')
# term values at t0
avals = [A0*t0, A1*(1-t0), B1*(1+kap1*t0), B4*(1+kap4*t0), B5]
logs = [mp.log(v) for v in avals]
# evaluate E3n(t0) with these logs
val = mp.mpf(0)
for key,poly in E3n.items():
    mval = poly[0]+poly[1]*t0+poly[2]*t0**2+poly[3]*t0**3
    lv = mp.mpf(1)
    for i,e in enumerate(key): lv *= logs[i]**e
    val += mval*lv
E3_engine = L*val
# direct d-series: |x_p'| at t0
y = mp.matrix([-h, halfL*(1-2*t0)])
yp = mp.matrix([mp.mpf(0), -2*halfL])
a_t = avals
ap_t = [A0, -A1, B1*kap1, B4*kap4, mp.mpf(0)]
n=3
def pmul(P,Q):
    R=[mp.mpf(0)]*(n+1)
    for i in range(n+1):
        for j in range(n+1-i): R[i+j]+=P[i]*Q[j]
    return R
def pdiv(P,Q):
    R=[mp.mpf(0)]*(n+1)
    for i in range(n+1): R[i]=(P[i]-sum(R[j]*Q[i-j] for j in range(i)))/Q[0]
    return R
def plog(P):
    Pp=[P[i]*i for i in range(1,n+1)]+[mp.mpf(0)]
    Q=pdiv(Pp,P); Ls=[mp.mpf(0)]
    for i in range(n): Ls.append(Q[i]/(i+1))
    return Ls
def pexp(P):
    E=[mp.mpf(1)]+[mp.mpf(0)]*n
    for k in range(1,n+1): E[k]=sum(P[j]*j*E[k-j] for j in range(1,k+1))/k
    return E
F=[mp.mpf(0)]*(n+1); S=[mp.mpf(0)]*(n+1)
for aj,lj in zip(a_t,logs):
    pk=mp.mpf(1)
    for k in range(n+1):
        if k: pk=pk*lj/k
        F[k]+=aj*pk
for aj,lj in zip(ap_t,logs):
    pk=mp.mpf(1)
    for k in range(n+1):
        if k: pk=pk*lj/k
        S[k]+=aj*pk
Phi=pdiv(S,F)
lg=plog(F)
inv=[mp.mpf(1),mp.mpf(-1),mp.mpf(1),mp.mpf(-1)]
m2lg=[-2*x for x in lg]
fac=pexp(pmul(m2lg,inv))
yp2=yp[0]**2+yp[1]**2; yyp=y[0]*yp[0]+y[1]*yp[1]; yy2=y[0]**2+y[1]**2
Phisq=pmul(Phi,Phi)
Q=[mp.mpf(0)]*(n+1); Q[0]=yp2
for k in range(1,n+1): Q[k]=-2*yyp*Phi[k]
for k in range(n+1): Q[k]+=yy2*Phisq[k]
G=pmul(fac,Q)
sq=[mp.mpf(0)]*(n+1); sq[0]=mp.sqrt(G[0])
for k in range(1,n+1): sq[k]=(G[k]-sum(sq[j]*sq[k-j] for j in range(1,k)))/(2*sq[0])
E3_direct = sq[3]
print("E3(0.3) engine:", mp.nstr(E3_engine,40))
print("E3(0.3) direct:", mp.nstr(E3_direct,40))
print("diff:", mp.nstr(abs(E3_engine-E3_direct),6))
with open('progress.log','a') as f:
    f.write(f"step: E3(0.3) check vs direct d-series: diff = {mp.nstr(abs(E3_engine-E3_direct),4)} (80 dps)\n")