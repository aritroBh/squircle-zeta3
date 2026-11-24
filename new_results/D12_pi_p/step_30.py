
import os
d = "<project>/research/08_theory/D12_pi_p"
os.makedirs(d, exist_ok=True)

files = {}

files['step_01.py'] = '''# step_01.py -- TASK 1: verify the integral formula and the duality pi_p = pi_q (60 digits)
import mpmath as mp
mp.mp.dps = 70
pi = mp.pi

def pi_raw(p):   # definition: pi_p = 2e int_0^1 (t^(1-p)+(1-t)^(1-p))^(1/p) dt, e=1/p
    p = mp.mpf(p); e = 1/p
    return 2*e*mp.quad(lambda t: (t**(1-p)+(1-t)**(1-p))**(1/p), [0, mp.mpf('0.5'), 1])

def pi_cosh(p):  # derived cosh form (proved in writeup): pi_p = 2h int_0^inf cosh(u/2)^(-1-h) cosh(u(1-h)/(2h))^h du
    p = mp.mpf(p); h = 1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])

print("pi_2 - pi =", mp.nstr(abs(pi_raw(2)-pi),3))
print("pi_1 =", pi_raw(1), " pi_inf ~ pi_1e6 =", mp.nstr(pi_cosh(mp.mpf('1e6')),15))
print("raw vs cosh:", [(p, mp.nstr(abs(pi_raw(p)-pi_cosh(p)),2)) for p in [2,3,mp.mpf('1.5'),5]])
print("duality table (p, q=conjugate, pi_p, |pi_p-pi_q|):")
for p in [mp.mpf('1.25'), mp.mpf('1.5'), 2, mp.mpf('2.5'), 3, 4, 5, 6, 10]:
    q = 1/(1-1/p); a,b = pi_cosh(p), pi_cosh(q)
    print(f"  p={mp.nstr(p,4)} q={mp.nstr(q,6)}: {mp.nstr(a,45)} |diff|={mp.nstr(abs(a-b),2)}")
print("min at p=2:", [(mp.nstr(p,5), mp.nstr(pi_cosh(p),20)) for p in [mp.mpf('1.9'),mp.mpf('1.95'),2,mp.mpf('2.05'),mp.mpf('2.1')]])
'''

files['step_02.py'] = '''# step_02.py -- TASK 2 setup: exact split pi_p = 4*2^(-1/p) + R(p), and J-integral evaluations
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi; z3 = mp.zeta(3)

def pi_cosh(p):
    p = mp.mpf(p); h = 1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])

def split(p):  # 4*2^(-h) + R(h), R via s-integral (Lemma 2 of writeup)
    p = mp.mpf(p); h = 1/p
    A = 4*mp.power(2,-h)
    w = lambda s: mp.exp(-s*h*h/(1-h))*(1+mp.exp(-s*h/(1-h)))**(-1-h)
    R = 4*h*h/(1-h)*mp.quad(lambda s: w(s)*((1+mp.exp(-s))**h-1), [0, mp.inf])
    return A, R

for p in [2,3,7,20]:
    A,R = split(p)
    print(f"p={p}: |(4*2^-h + R) - pi_p| = {mp.nstr(abs(A+R-pi_cosh(p)),2)}")

# hypergeometric remark: g(h) = 4h int_0^1 v^(h-1)(1+v)^(-1-h) dv = 4*2F1(1+h,h;1+h;-1) = 4*2^(-h)
g = lambda h: 4*h*mp.quad(lambda v: v**(h-1)*(1+v)**(-1-h), [0,1])
for h in [mp.mpf('0.3'), mp.mpf('0.77')]:
    print(f"g({h}) - 4*2^(-h) = {mp.nstr(abs(g(h)-4*mp.power(2,-h)),2)}")

# J(m,k) = int_0^inf s^k log(1+e^-s)^m ds : key evaluations
def J(m,k): return mp.quad(lambda s: s**k*mp.log(1+mp.e**(-s))**m, [0, mp.inf])
print("J(1,0)=pi^2/12:", mp.nstr(abs(J(1,0)-pi**2/12),2))
print("J(1,1)=3zeta(3)/4:", mp.nstr(abs(J(1,1)-3*z3/4),2))
print("J(2,0)=zeta(3)/4:", mp.nstr(abs(J(2,0)-z3/4),2))
'''

files['step_03.py'] = '''# step_03.py -- TASK 2: closed forms N0..N3 (PROVED in writeup) and numerical confirmation
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi; L2 = mp.log(2); z3 = mp.zeta(3)
N0 = -L2
N1 = L2**2/2
N2 = pi**2/24 - L2**3/6
N3 = L2**4/24 + pi**2/24*(1-L2) + z3/4

def pi_cosh(p):
    p = mp.mpf(p); h = 1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])

print("N0..N3 =", [mp.nstr(x,40) for x in [N0,N1,N2,N3]])
for p in [40, 80, 160]:
    h = 1/mp.mpf(p); f = pi_cosh(p)/4
    S3 = 1 + N0*h + N1*h**2 + N2*h**3
    S4 = S3 + N3*h**4
    print(f"p={p}: f-S3={mp.nstr(f-S3,3)} (O(h^5), h^5={mp.nstr(h**5,3)})   f-S4={mp.nstr(f-S4,3)}")
'''

files['step_04.py'] = '''# step_04.py -- TASK 3: exact series arithmetic for [h^n] Psi (rational coeffs in s, L=log(1+e^-s), log2)
from fractions import Fraction
import math, pickle
import sympy as sp
KMAX = 44
def addC(A,B):
    C=dict(A)
    for k,v in B.items():
        C[k]=C.get(k,Fraction(0))+v
        if C[k]==0: del C[k]
    return C
def scaleC(A,c): return {k:v*c for k,v in A.items()}
def mulC(A,B):
    C={}
    for (i1,j1,k1),v1 in A.items():
        for (i2,j2,k2),v2 in B.items():
            key=(i1+i2,j1+j2,k1+k2); C[key]=C.get(key,Fraction(0))+v1*v2
    return {k:v for k,v in C.items() if v!=0}
def serMul(A,B,K=KMAX):
    C=[{} for _ in range(K+1)]
    for i in range(len(A)):
        if not A[i]: continue
        for j in range(len(B)):
            if not B[j]: continue
            if i+j>K: break
            C[i+j]=addC(C[i+j],mulC(A[i],B[j]))
    return C
def serExp(A,K=KMAX):
    R=[{} for _ in range(K+1)]; R[0]={(0,0,0):Fraction(1)}
    for n in range(1,K+1):
        acc={}
        for k in range(1,n+1):
            if A[k]: acc=addC(acc,mulC(scaleC(A[k],Fraction(k)),R[n-k]))
        R[n]=scaleC(acc,Fraction(1,n))
    return R
y=sp.symbols('y')
lc=sp.log(sp.cosh(y)).series(y,0,KMAX+2).removeO()
c=[None]*(KMAX+1); c[1]=Fraction(-1,2)
for r in range(2,KMAX+1): c[r]=Fraction(sp.expand(lc).coeff(y,r))/Fraction(2)**r
T=[{} for _ in range(KMAX+1)]; T[0]={(0,0,1):Fraction(1)}   # log(1+e^{-sa}) = log2 + sum c_r s^r a^r
for r in range(1,KMAX+1):
    if c[r]==0: continue
    for t in range(0,KMAX+1-r):
        T[r+t]=addC(T[r+t],{(r,0,0):Fraction(math.comb(r+t-1,t))*c[r]})
Q=[{} for _ in range(KMAX+1)]
for n in range(KMAX+1):
    Q[n]=addC(Q[n],scaleC(T[n],Fraction(-1)))
    if n+1<=KMAX: Q[n+1]=addC(Q[n+1],scaleC(T[n],Fraction(-1)))
psi2=[scaleC(co,Fraction(1,2)) for co in serExp([{}]+Q[1:])]   # (1+e^{-sa})^{-1-h}, exp(Q0)=1/2
B=[{} for _ in range(KMAX+1)]
for m in range(2,KMAX+1): B[m]={(1,0,0):Fraction(-1)}
psi1=serExp(B)                                                  # exp(-s h^2/(1-h))
psi3=[{} for _ in range(KMAX+1)]
for m in range(1,KMAX+1): psi3[m]={(0,m,0):Fraction(1,math.factorial(m))}  # exp(hL)-1
PSI=serMul(serMul(psi1,psi2),psi3)
pickle.dump({'PSI':PSI,'KMAX':KMAX},open('/tmp/psi.pkl','wb'))
print("low orders:", PSI[1], PSI[2])
print("saved PSI to /tmp/psi.pkl")
'''

files['step_05.py'] = '''# step_05.py -- TASK 3: assemble N_k (k<=43) from PSI via J-integrals at 110 digits; cross-check by direct fit
import mpmath as mp, pickle
mp.mp.dps = 110
PSI = pickle.load(open('/tmp/psi.pkl','rb'))['PSI']; KMAX=44
pairs=set()
for co in PSI:
    for (i,j,k) in co: pairs.add((j,i))
J={}
for (j,i) in sorted(pairs):
    J[(j,i)] = mp.quad(lambda s: s**i*mp.log(1+mp.e**(-s))**j, [0, mp.inf])
log2=mp.log(2)
INT=[mp.mpf(0)]*(KMAX+1)
for m in range(KMAX+1):
    INT[m] = sum(mp.mpf(v.numerator)/mp.mpf(v.denominator)*J[(j,i)]*log2**k
                 for (i,j,k),v in PSI[m].items())
N=[mp.mpf(0)]*(KMAX+2)
for k in range(0,KMAX+1):
    N[k] = (-log2)**(k+1)/mp.factorial(k+1) + sum(INT[0:k])
for k in range(12): print(f"N_{k:2d} = {mp.nstr(N[k],40)}")
pickle.dump({'N':N,'INT':INT,'J':J}, open('/tmp/Nvals.pkl','wb'))
# cross-check first 11 coefficients by direct fit of pi_{1/h}/4 at 150 digits
mp.mp.dps=150
def pi_h(h):
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f,[0,mp.inf])/4
ps=[60+3*i for i in range(14)]; fs=[pi_h(1/mp.mpf(p)) for p in ps]
A=mp.matrix(14,11); b=mp.matrix(14,1)
for a,p in enumerate(ps):
    h=1/mp.mpf(p)
    for k in range(11): A[a,k]=h**(k+1)
    b[a]=fs[a]-1
coef=mp.lu_solve(A,b)
mp.mp.dps=40
print("direct fit vs series arithmetic:", [mp.nstr(abs(coef[k]-N[k]),2) for k in range(11)])
'''

files['step_06.py'] = '''# step_06.py -- TASK 2/4: PSLQ identifications (95 digits) and sqrt2-absence test
import mpmath as mp, pickle
from math import gcd
mp.mp.dps = 95
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3); z5=mp.zeta(5)
Li4=mp.polylog(4,mp.mpf('0.5')); Li5=mp.polylog(5,mp.mpf('0.5'))
def identify(tag,t,b,nm):
    rel=mp.pslq([t]+b,maxcoeff=10**10,maxsteps=10000)
    g=0
    for x in rel: g=gcd(g,x)
    rel=[x//g for x in rel]
    val=-sum(rel[i+1]*b[i] for i in range(len(b)))/rel[0]
    print(tag, rel, "residual=", mp.nstr(abs(val-t),2))
identify('N2',N[2],[L2**3,pi**2],[])
identify('N3',N[3],[L2**4,pi**2,pi**2*L2,z3],[])
identify('N4',N[4],[L2**5,pi**2,pi**2*L2,pi**2*L2**2,z3,z3*L2,pi**4],[])
# J-evaluations at weights 4 and 5 (Euler sums)
def J(m,k): return mp.quad(lambda s: s**k*mp.log(1+mp.e**(-s))**m,[0,mp.inf])
W4=[L2**4,pi**2*L2**2,pi**4,z3*L2,Li4]
identify('J30',J(3,0),W4,[]); identify('J21',J(2,1),W4,[])
W5=[L2**5,pi**2*L2**3,pi**4*L2,z3*L2**2,z3*pi**2,z5,Li4*L2,Li5]
identify('J40',J(4,0),W5,[]); identify('J31',J(3,1),W5,[]); identify('J22',J(2,2),W5,[])
# N5 closed form check
N5form=(L2**6/720 + pi**2/24*(1-L2+L2**2/2-L2**3/6) + z3/4*(1-L2+L2**2/2)
        + pi**4/2880*(1-L2) - 7*pi**4/1920 - z3*pi**2/48 + z5/8)
print("N5 form residual:", mp.nstr(abs(N5form-N[5]),2))
# sqrt2 absence: give PSLQ BOTH alphabets
sq2=mp.sqrt(2); las=mp.log(1+sq2); li2m=mp.polylog(2,sq2-1)
bE=[sq2,las,sq2*las,li2m,las**2,L2**4,pi**2,pi**2*L2,z3]
rel=mp.pslq([N[3]]+bE,maxcoeff=10**9,maxsteps=10000)
print("N3 relation in joint basis (first 5 slots = Euclidean constants):", rel)
'''

files['step_07.py'] = '''# step_07.py -- TASK 3: growth analysis, radius, singularity asymptotics, optimal truncation
import mpmath as mp, pickle
mp.mp.dps = 80
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']; INT=pickle.load(open('/tmp/Nvals.pkl','rb'))['INT']
pi=mp.pi
print("k : N_k/prev  |N_k|^{1/k}  k/(e*pi)  sqrt(|N_{k+2}/N_k|)*pi/(k+1)  odd-norm N_k*pi^{k-1}/(2(k-2)!)")
for k in range(10,41):
    rt=abs(N[k])**(1/k); q=mp.sqrt(abs(N[k+2]/N[k]))*pi/(k+1)
    norm = N[k]*pi**(k-1)/(2*mp.factorial(k-2)) if k%2==1 else mp.nan
    print(f"{k:2d} {mp.nstr(N[k+1]/N[k],8):>10} {mp.nstr(rt,8)} {mp.nstr(k/(mp.e*pi),8)} {mp.nstr(q,8)}  {mp.nstr(norm,8) if k%2==1 else ''}")
# integrand singularity asymptotics: [h^m]psi2 ~ -(2/s) Re(1+i s/pi)^{m+1}
def psi2(s,h): return (1+mp.e**(-s*h/(1-h)))**(-1-h)
def coeff_psi2(s,m):
    rho=0.8*pi/mp.sqrt(s**2+pi**2)
    return mp.quad(lambda t: psi2(s,rho*mp.e**(1j*t))/(rho*mp.e**(1j*t))**m,[0,2*pi])/(2*pi)
for s in [10.0,30.0]:
    for m in [9,10]:
        cm=coeff_psi2(s,m); asym=-(1/s)*2*mp.re((1+1j*s/pi)**(m+1))
        print(f"s={s} m={m}: [h^m]psi2={mp.nstr(cm,10)} asym={mp.nstr(asym,10)} ratio={mp.nstr(cm/asym,6)}")
# optimal truncation at p=10
mp.mp.dps=110
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']
def pi_cosh(p):
    p=mp.mpf(p); h=1/p
    f=lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f,[0,mp.inf])
true=pi_cosh(10)/4; h=mp.mpf('0.1'); s=mp.mpf(1); best=(0,mp.mpf(1))
for k in range(0,42):
    s+=N[k]*h**(k+1); err=abs(s-true)
    if err<best[1]: best=(k,err)
print("optimal truncation p=10: best k=",best[0]," err=",mp.nstr(best[1],3),
      " pi/h=",mp.nstr(pi/h,5)," e^(-pi/h)=",mp.nstr(mp.e**(-pi/h),3))
'''

files['step_08.py'] = '''# step_08.py -- TASK 3 (p=1 expansion) and TASK 4 (Euclidean comparison)
import mpmath as mp, pickle
mp.mp.dps = 110
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']; KMAX=44
f=[mp.mpf(1)]+[N[k] for k in range(KMAX+1)]
# (a) task variable eta = 1/q - 1 = -h :  pi_q/4 = 1 + sum_{k>=1} (-1)^k N_{k-1} eta^k
print("task-3 series in eta=1/q-1: coefficients (-1)^k N_{k-1}:")
for k in range(1,13): print(f"  {(-1)**k}*N_{k-1} = {mp.nstr((-1)**k*N[k-1],35)}")
# (b) alternative variable eta=q-1: h=eta/(1+eta) -> M_n
M=[mp.mpf(0)]*(KMAX+2); M[0]=mp.mpf(1)
for n in range(1,KMAX+2):
    M[n]=sum(f[j]*(-1)**(n-j)*mp.binomial(n-1,j-1) for j in range(1,n+1))
print("M_n (pi_{1+eta}/4) for n=1..15:", [mp.nstr(M[n],12) for n in range(1,16)])
print("growth |M_n|^(1/n):", [mp.nstr(abs(M[n])**(1/n),8) for n in [10,15,20,25,30,35,39]])
# verify against direct evaluation (optimal truncation)
def pi_cosh(p):
    p=mp.mpf(p); h=1/p
    g=lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(g,[0,mp.inf])
for eta in [mp.mpf('0.02'),mp.mpf('0.05'),mp.mpf('0.1')]:
    true=pi_cosh(1+eta)/4; s=mp.mpf(0); best=(0,mp.mpf(1))
    for n in range(0,42):
        s+=M[n]*eta**n; err=abs(s-true)
        if err<best[1]: best=(n,err)
    print(f"eta={eta}: best trunc n={best[0]} err={mp.nstr(best[1],3)} (pi/h with h=eta/(1+eta): {mp.nstr(mp.pi*(1+eta)/eta,6)})")
# TASK 4: Euclidean constants
print("\\nEuclidean M0 = log2 - sqrt2*log(1+sqrt2) =", mp.nstr(mp.log(2)-mp.sqrt(2)*mp.log(1+mp.sqrt(2)),40))
print("our N0 = -log2 =", mp.nstr(-mp.log(2),40))
'''

for name, content in files.items():
    with open(os.path.join(d,name),'w') as fh:
        fh.write(content)
    print("wrote", name)
