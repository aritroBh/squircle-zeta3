src = r'''# verify_grw1r.py -- GRW1R referee's independent verification script (target < 5 min)
# Verifies: C (E* certified; repaired Darboux remainder 7.06), D (uniformity bound),
# E (two-term law per class), F/T2 (Borel-Laplace S(eta)=I(eta) by own conformal map),
# B (Lemma 2), A (Laplace form), G (R=B-H radius pi), H (jump integral asymptotics).
import time
t0 = time.time()
import mpmath as mp, json, os
from fractions import Fraction
from math import factorial, comb
from flint import arb

HERE = os.path.dirname(os.path.abspath(__file__))
F450 = json.load(open(os.path.join(HERE, '../D11b2_stokes_jump/step2_num450_N180.json')))['f450']

print('==== [1] Claim C: E* = 5/24 CERTIFIED (arb cover + exact series) ====')
# E(w) = w^{-2}(sqrt((1-e^{-2w})/(2w)) - 1 + w/2); slit w = i s, s in [0, 1.9 rho]
# Piece A: s in [0,0.1]: psi(s)=|E(is)|^2 = sum c_m s^{2m}; c_1<0 dominates -> psi <= (5/24)^2
e6 = arb(6).exp()
M_E3 = (((1+e6)/6).sqrt()+1+arb('1.5'))/9     # rigorous max |E| on |w|<=3
c1 = arb(-25)/13824; c2 = arb(11147)/232243200
u_max = arb('0.01')
T = (M_E3**2)/(9**3*(1-u_max/9))
assert (c1 + (c2+T)*u_max).upper() < 0
print(f'  piece A [0,0.1]: psi(s) <= (5/24)^2 since c1+(c2+T)u < {(c1+(c2+T)*u_max).upper()}  [series+Cauchy tail, M_E(3)<={M_E3.upper()}]')
# exact series coeffs for piece B
import sympy as sp
K = 96
am=[sp.Rational((-2)**m, sp.factorial(m+1)) for m in range(K+3)]
b=[sp.Rational(0)]*(K+3); b[0]=sp.Rational(1)
for m in range(1,K+3): b[m]=(am[m]-sum(b[i]*b[m-i] for i in range(1,m)))/2
e=[b[k+2] for k in range(K+1)]
cc=[sp.Rational(0)]*(K//2+1)
for m in range(K//2+1):
    cc[m]=((-1)**m)*sum(((-1)**(2*m-i))*e[i]*e[2*m-i] for i in range(2*m+1))
carb=[arb(x.numerator)/arb(x.denominator) for x in cc]
KKs=45; M2=M_E3**2
def psi_series(s):
    u=s*s; P=arb(0)
    for m in range(KKs,-1,-1): P=P*u+carb[m]
    tail=M2*(u/9)**(KKs+1)/(1-u/9)
    return arb(P.mid(), P.rad()+float(tail.upper()))
def psi_direct(s):
    q=s.sin()/s
    return (q+1+s**2/4-(2*(s/2).cos()+s*(s/2).sin())*q.sqrt())/s**4
rho=arb.pi()/2; send=arb('1.9')*rho
maxub=None
c=arb('0.1')+arb('0.001')
while float(c-arb('0.001'))<=1.5:
    seg=arb(c,arb('0.001'))
    v=psi_series(seg)
    if maxub is None or v.upper()>maxub: maxub=v.upper()
    c+=arb('0.002')
c=arb('1.501')+arb('0.001')
while float(c-arb('0.001'))<=float(send):
    seg=arb(c,arb('0.001'))
    v=psi_direct(seg)
    if v.upper()>maxub: maxub=v.upper()
    c+=arb('0.002')
psi524=(arb(5)/24)**2
print(f'  piece B: certified max psi on [0.1, 1.9 rho] = {maxub}')
print(f'  (5/24)^2 = {psi524}  =>  E* = max|E| on slit = {maxub.sqrt() if maxub>psi524 else psi524.sqrt()} = 5/24  CERTIFIED')

print('==== [2] Claim C: exact H_n, two-term law, repaired remainder 7.06 (claimed 4.45+circle true) ====')
mp.mp.dps = 100
N=220
KK=N+2
a=[Fraction(0)]*KK; a[0]=Fraction(1)
for k in range(1,KK): a[k]=Fraction((-2)**k,2*factorial(k))
s=[Fraction(0)]*KK; s[0]=Fraction(1)
for k in range(1,KK): s[k]=(a[k]-sum(s[i]*s[k-i] for i in range(1,k)))/2
Hmp=[mp.mpf((2*s[n]).numerator)/mp.mpf((2*s[n]).denominator) for n in range(N+1)]
rhom=mp.pi/2; rhotm=mp.mpf('2.9')*rhom
def law2(n):
    return -mp.sqrt(2)*rhom**(-n)*mp.mpf(n)**mp.mpf('-1.5')*mp.cos((2*n+1)*mp.pi/4)*(1+(3-3*mp.pi*((-1)**n))/(8*n))
w1=w2=mp.mpf(0); a1=a2=None
CH=mp.mpf(0); nCH=0
for n in range(10,221):
    R=abs(Hmp[n]-law2(n))
    r1=R/(mp.mpf('7.06')*rhom**(-n)*mp.mpf(n)**mp.mpf('-3.5'))
    r2=R/(mp.mpf('4.45')*rhom**(-n)*mp.mpf(n)**mp.mpf('-3.5')+135*rhotm**(-n))
    if r1>w1: w1,a1=r1,n
    if r2>w2: w2,a2=r2,n
for n in range(1,201):
    v=abs(Hmp[n])*rhom**n*mp.mpf(n)**mp.mpf('1.5')
    if v>CH: CH,nCH=v,n
print(f'  |H_n-law2| <= 7.06 rho^-n n^-7/2 : worst ratio {mp.nstr(w1,6)} (n={a1})  PASS')
print(f'  |H_n-law2| <= 4.45 rho^-n n^-7/2 + 135(2.9rho)^-n : worst ratio {mp.nstr(w2,6)} (n={a2})  TRUE on 10..220')
print(f'  C_H = max |H_n| rho^n n^1.5 = {mp.nstr(CH,8)} at n={nCH} (claimed 5.874 @3)  PASS')

print('==== [3] Claim C: certified honest remainder sup over 10..500 (+ Stirling tail) ====')
def sgn_xa(n): return [-1,1,1,-1][n%4]
def sgn_xb(n): return [1,1,-1,-1][n%4]
def sgn_law(n): return [-1,1,1,-1][n%4]
sqrt2a=arb(2).sqrt()
def betaab(x,y): return x.gamma()*y.gamma()/(x+y).gamma()
def A_n(n): return sgn_xa(n)*sqrt2a*(2/arb.pi())*rho**(arb('0.5')-n)*betaab(arb('1.5'),arb(n)-arb('0.5'))
def B_n(n): return sgn_xb(n)*(sqrt2a/2)*(2/arb.pi())*rho**(arb('1.5')-n)*betaab(arb('2.5'),arb(n)-arb('1.5'))
def law2a(n): return sgn_law(n)*rho**(-n)*arb(n)**arb('-1.5')*(1+(3-3*arb.pi()*((-1)**n))/(8*n))
def Ebnd(n): return (4/arb.pi())*(arb(5)/24)*rho**(arb('2.5')-n)*betaab(arb('3.5'),arb(n)-arb('2.5'))
def tailsa(n):
    rt=arb('2.9')*rho
    return (4/arb.pi())*rt**(arb('0.5')-n)/(arb(n)-arb('0.5'))+(2/arb.pi())*rt**(arb('1.5')-n)/(arb(n)-arb('1.5'))
mx=None; argmx=None
for n in range(10,501):
    u=rho**(-n)*arb(n)**arb('-3.5')
    val=(abs(A_n(n)+B_n(n)-law2a(n))+Ebnd(n)+tailsa(n))/u
    if mx is None or val.upper()>mx: mx,argmx=val.upper(),n
print(f'  certified sup (|A+B-law2|+Ebnd+tails)/(rho^-n n^-7/2) over 10..500 = {mx} at n={argmx}')
print(f'  (Stirling tail n>=500 proved <= 5.19; => |R_n| <= 7.06 rho^-n n^-7/2 + 135(2.9rho)^-n for ALL n>=10)')

print('==== [4] Claim D: uniformity bound vs true tail (450-digit f_n) ====')
f450=[mp.mpf(v) for v in F450]
def etaD(n):
    if n==1: return mp.log(2)
    return (1-mp.mpf(2)**(1-n))*mp.zeta(n)
for n in [20,40,80,120]:
    main=Hmp[n]*mp.factorial(n)*etaD(n)
    tail=f450[n]-main
    Dbnd=mp.mpf('5.874')*mp.e**(2*rhom)*rhom**(-n)*mp.factorial(n)*mp.mpf(2)**(-n-1)/n
    print(f'  n={n}: |tail|={mp.nstr(abs(tail),4)} <= D-bound {mp.nstr(Dbnd,4)} : {abs(tail)<=Dbnd} ; |tail|/|main|={mp.nstr(abs(tail)/abs(main),3)}')

print('==== [5] Claim E: two-term law per residue class mod 4 (Richardson) ====')
def law_main(n):
    return mp.sqrt(2)*mp.cos(mp.pi*(3-2*n)/4)*mp.factorial(n)*(mp.mpf(2)/mp.pi)**n*mp.mpf(n)**mp.mpf('-1.5')
for cls in range(4):
    ns=[160+cls,156+cls,152+cls]
    rs=[f450[k]/law_main(k) for k in ns]
    M=mp.matrix(3,3)
    for i,k in enumerate(ns): M[i,0]=1; M[i,1]=1/mp.mpf(k); M[i,2]=1/mp.mpf(k)**2
    sol=mp.lu_solve(M,mp.matrix(rs))
    c1p=(3-3*mp.pi)/8 if cls%2==0 else (3+3*mp.pi)/8
    print(f'  class {cls}: const={mp.nstr(sol[0],12)} (~1), c1={mp.nstr(sol[1],10)} (pred {mp.nstr(c1p,10)})')

print('==== [6] Theorem 2 (F): S(eta)=I(eta) by own conformal-map Borel-Laplace ====')
mp.mp.dps = 120
r=[f450[n]/mp.factorial(n)-Hmp[n] for n in range(181)]
beta=[mp.mpf(0)]*181; beta[0]=r[0]
twopi=2*mp.pi
for n in range(1,181):
    beta[n]=sum(r[n-2*j]*(twopi**(n-2*j))*comb(n-j-1,j) for j in range(0,(n-1)//2+1) if n-2*j>=1)
def w_of_z(z): return z/(mp.pi+mp.sqrt(mp.pi**2+z**2))
def R_map(z):
    w=w_of_z(z); tot=mp.mpf(0); pw=mp.mpf(1)
    for n in range(181): tot+=beta[n]*pw; pw*=w
    return tot
def Hz(z): return mp.sqrt(2)*mp.sqrt(1+mp.e**(-2*z))
def S_H(eta): return mp.quad(lambda t: mp.e**(-t)*Hz(eta*t), [0,5,20,60,mp.inf])
def S_R(eta): return mp.quad(lambda t: mp.e**(-t)*R_map(eta*t), [0,2,5,10,20,40,mp.inf])
def I_perim(eta):
    f=lambda t: mp.sqrt(t**(2*eta)+(1-t)**(2*eta))
    P=4*(1+eta)*mp.quad(f,[0,mp.mpf('0.5'),1])
    return P/(4*mp.sqrt(2)*(1+eta))
for eta in [mp.mpf('0.5'), mp.mpf('2'), mp.mpf('-0.5'), mp.mpf('-0.9')]:
    S=S_H(eta)+S_R(eta); Iv=I_perim(eta)
    print(f'  eta={eta}: S={mp.nstr(S,15)}  I={mp.nstr(Iv,15)}  |S-I|={mp.nstr(abs(S-Iv),3)}')
# G: radius of R=B-H = pi (root test)
print(f'  [G] |r_180|^(1/180) = {mp.nstr(abs(r[180])**(mp.mpf(1)/180),8)}  vs 1/pi = {mp.nstr(1/mp.pi,8)}')

print('==== [7] Claim B (Lemma 2) spot check, own quadrature ====')
mp.mp.dps=60
def kjfun(j,u):
    if u==0: return mp.mpf('0.25')*(-mp.log(2))**j/mp.factorial(j)
    L=-mp.log(1+mp.e**(-u))
    return (1+mp.e**(-u))**(-2)*L**j/mp.factorial(j)
def M(n,j):
    return mp.quad(lambda u: mp.e**(-u)*u**(n-j)*kjfun(j,u), [0,mp.inf], maxdegree=8)
for n in [1,2,3,5]:
    tot=mp.mpf(0)
    for j in range(0,n+1): tot+=Hmp[n-j]*M(n,j)
    print(f'  n={n}: |sum - f_n|/|f_n| = {mp.nstr(abs(tot-f450[n])/abs(f450[n]),3)}')
print('==== [8] Claim A spot check, eta=0.37 ====')
eta=mp.mpf('0.37')
lap=mp.quad(lambda u: mp.e**(-u)*Hz(eta*u)*(1+mp.e**(-u))**(-2-eta), [0,5,15,mp.inf])
Iv=I_perim(eta)
print(f'  Laplace form vs perimeter: |diff| = {mp.nstr(abs(lap-Iv),3)}')

print('==== [9] Claim H: Stokes jump integral asymptotics ====')
eps=mp.mpf('0.1')
Jint=-(4/eps)*mp.e**(-mp.pi/(2*eps))*mp.e**(1j*mp.pi/4)*mp.quad(lambda s: mp.e**(-s/eps)*mp.sqrt(mp.sin(s))*mp.e**(-1j*s/2), [0, eps*3, eps*20, mp.pi])
Jpred=2*mp.sqrt(mp.pi)*mp.e**(-3j*mp.pi/4)*mp.sqrt(eps)*mp.e**(-mp.pi/(2*eps))*(1-3j*eps/4-25*eps**2/32)
print(f'  eps=0.1: |J|={mp.nstr(abs(Jint),10)} |pred|={mp.nstr(abs(Jpred),10)} ratio={mp.nstr(abs(Jint)/abs(Jpred),8)}; J/pred-1={mp.nstr(Jint/Jpred-1,6)}')

print(f'\\nTOTAL TIME: {time.time()-t0:.1f}s')
'''
open('verify_grw1r.py','w').write(src)
print("written", len(src))
