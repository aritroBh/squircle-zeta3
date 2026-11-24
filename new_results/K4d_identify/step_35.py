import mpmath as mp, json
from mpmath import mpf
mp.mp.dps = 70
WORK='research/08_theory/K4d_identify'
SA = [mpf(s) for s in json.load(open(WORK+'/modelA_h96.json'))['Sdk']]
S0,S1,S2,S3,S4 = SA
# c_A = 2^(-6/5) 3^(13/5); L = ln c_A
L = (mpf(-6)/5)*mp.log(2) + (mpf(13)/5)*mp.log(3)
# T-series: S(1+d) = c_A^{d/(1+d)} T(d); factor f(d)=c_A^{d/(1+d)}
# f = exp(L*u), u=d/(1+d)=d-d^2+d^3-d^4+...
# f coeffs: f0=1, f1=L, f2=L^2/2 - L, f3=L - L^2 + L^3/6, f4=...
u = [mpf(0)]*5
u[1]=1; u[2]=-1; u[3]=1; u[4]=-1   # u = d/(1+d) coeffs
# exp(L*u): e = sum (Lu)^k/k!
import itertools
def sermul(a,b,K):
    o=[mpf(0)]*(K+1)
    for i in range(K+1):
        for j in range(i+1): o[i]+=a[j]*b[i-j]
    return o
Lu=[L*x for x in u]
e=[mpf(0)]*5; e[0]=1
pw=[mpf(0)]*5; pw[0]=1
import math
for k in range(1,5):
    pw=sermul(pw,Lu,4)
    for i in range(4,-1,-1):
        e[i]+=pw[i]/mp.factorial(k)
# f = e; S = f * T  =>  T = S / f  (series division)
def serinv(a,K):
    o=[mpf(0)]*(K+1); o[0]=1/a[0]
    for n in range(1,K+1):
        o[n]=-sum(a[j]*o[n-j] for j in range(1,n+1))/a[0]
    return o
finv=serinv(e,4)
T=sermul(list(SA),finv,4)
print("factor coeffs f:", [mp.nstr(x,12) for x in e])
print("T0 =", mp.nstr(T[0],50))
print("T1 =", mp.nstr(T[1],50), "  T1/S0=",mp.nstr(T[1]/S0,30), " (13/15=%.15f)"%(13/15))
print("T2 =", mp.nstr(T[2],50))
print("T2/S0 =", mp.nstr(T[2]/S0,55))
print("T3 =", mp.nstr(T[3],50))
print("T3/S0 =", mp.nstr(T[3]/S0,55))
print("T4 =", mp.nstr(T[4],40))
print("T4/S0 =", mp.nstr(T[4]/S0,45))
