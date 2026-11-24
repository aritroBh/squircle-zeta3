# verify_pwr1.py -- PWR1 referee verification battery (self-contained; runtime < 5 min).
# Verifies: (A) thm:pw exact mechanism + numerics; (B) prop:degree machinery + numerics.
# Artifacts used: pwr1_engine.py, pwr1_rhombus.py (referee's engines, written during review),
#                 ../D16_pure_weight/tau_series_2d.json, ../K3_fourier_all_orders/e_rm.pkl, Qpolys.pkl
import mpmath as mp, sympy as sp, json, pickle, time
from fractions import Fraction
t00 = time.time()
mp.mp.dps = 70
PASS = []
def report(name, ok, detail=""):
    PASS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")

print("="*80); print("PART A: thm:pw (pure weight)"); print("="*80)
import pwr1_engine as E
E.NMAX = 8
tot2, pure2, B2 = E.total_log(2)
gam_ok = all(sp.expand(tot2.coeff(k)).coeff(E.GAM,1)==0 for k in range(9))
w_ok = all(E.weight_check(pure2.coeff(k), k) for k in range(9))
report("A1 gamma cancels exactly, n=2, all k<=8", gam_ok)
report("A2 weight-k homogeneity of w_k(2), all k<=8", w_ok)
report("A3 Rat_2 = 1 (logRat vanishes)", E.logRat(2).d == {})
u = [sp.expand(pure2.coeff(k)).coeff(E.Z[k],1) for k in range(2,9)]
report("A4 single-zeta u_k(2) = 1/4,1/4,1/32,-3/16,91/192,129/64,-5039/1024",
       u == [sp.Rational(s) for s in ['1/4','1/4','1/32','-3/16','91/192','129/64','-5039/1024']], str(u))
# exact vs D16 json
js = json.load(open('../D16_pure_weight/tau_series_2d.json'))
def mono(expr):
    out = {}
    for t in sp.expand(expr).as_ordered_terms() if expr != 0 else []:
        c, zp = t.as_coeff_Mul(); key = []
        for f in sp.Mul.make_args(zp):
            if f == 1: continue
            b, ee = f.as_base_exp()
            zi = [k for k,v in E.Z.items() if v==b][0]; key.extend([zi]*int(ee))
        out[' '.join(map(str,sorted(key)))] = out.get(' '.join(map(str,sorted(key))),0)+c
    return {k: sp.nsimplify(v) for k,v in out.items()}
P = E.ser_exp(pure2).scale(4)
okP = all(mono(P.coeff(k)) == {kk: sp.Rational(vv) for kk,vv in js['P_over_sqrt2'].get(str(k),{}).items()} for k in range(9))
okL = all(mono(pure2.coeff(k)) == {kk: sp.Rational(vv) for kk,vv in js['logP'].get(str(k),{}).items()} for k in range(9))
report("A5 P/sqrt2 + logP exact vs D16 tau_series_2d.json through tau^8", okP and okL)
# n=3 spot (fast): w2, w3, Rat_3
tot3, pure3, _ = E.total_log(3)
report("A6 w_2(3) = -2 Z2/3, w_3(3) = -16 Z3/9 (H1 T5 L2,L3)",
       pure3.coeff(2) == -sp.Rational(2,3)*E.Z[2] and pure3.coeff(3) == -sp.Rational(16,9)*E.Z[3])
lR3 = E.logRat(3)
report("A7 logRat_3 = log(1-tau)-log(1-2 tau)  [Rat_3 = 1/(2-p)]",
       all(lR3.coeff(k) == sp.Rational(2**k-1, k) for k in range(1,9)))
gam3 = all(sp.expand(tot3.coeff(k)).coeff(E.GAM,1)==0 for k in range(9))
w3ok = all(E.weight_check(pure3.coeff(k), k) for k in range(9))
report("A8 n=3: gamma-free, pure weight through tau^8", gam3 and w3ok)
print(f"   [engine part: {time.time()-t00:.0f}s]")

# ---- numeric remainder tests (fresh quadrature) ----
w = json.load(open('w_exact_n2_tau10.json'))
znum = {E.Z[k]: mp.zeta(k) for k in range(2,12)}
W = {int(k): mp.mpf(str(sp.sympify(v).subs(znum).evalf(50))) for k,v in w.items()}
def logP_quad(tau):
    tau = mp.mpf(tau); q = 1/tau
    h = lambda th: (mp.cos(th)**q + mp.sin(th)**q)**tau
    return mp.log(4*mp.quad(h, [0, mp.pi/4, mp.pi/2])/(4*mp.sqrt(2)))
tau = mp.mpf('0.002'); v = logP_quad(tau)
R = v - mp.fsum(W[k]*tau**k for k in range(2,9)); ratio = R/(W[9]*tau**9)
R2 = v - mp.fsum(W[k]*tau**k for k in range(2,10)); ratio2 = R2/(W[10]*tau**10)
print(f"   n=2 remainder test at tau=0.002: R/(w9 tau^9) = {mp.nstr(ratio,10)} (->1), (R-w9 t9)/(w10 t10) = {mp.nstr(ratio2,10)} (->1)")
report("A9 quadrature remainder test (expect ~0.988, ~1.011)", abs(ratio-1) < 0.02 and abs(ratio2-1) < 0.02)

# ---- rhombus ----
import pwr1_rhombus as RH
s = RH.logP_rhombus()
c2, c3 = sp.expand(s).coeff(RH.TAU,2), sp.expand(s).coeff(RH.TAU,3)
report("A10 rhombus w2 = (1-q^2) Z2/4, w3 = (1-q^2) Z3/4",
       sp.expand(c2 - (1-RH.q**2)*RH.Z[2]/4)==0 and sp.expand(c3 - (1-RH.q**2)*RH.Z[3]/4)==0)
t = sp.Symbol('t')
pur = True
for k in range(1,7):
    c = sp.expand(s).coeff(RH.TAU, k)
    sub = sp.expand(c.subs({RH.Z[i]: RH.Z[i]*t**i for i in RH.Z}))
    pur &= all(m[0]==k for m in sp.Poly(sub,t).monoms()) if sub!=0 else True
    pur &= all(m[0]%2==0 for m in sp.Poly(c,RH.q).monoms()) if c!=0 else True
report("A11 rhombus: pure weight + coefficients in Q[q^2] through tau^6 (exact)", pur)

print("="*80); print("PART B: prop:degree (degree of Fourier coefficients)"); print("="*80)
# e table rebuild (exact, fast)
ORD = 14
def smul(a,b,n=ORD):
    out={}
    for i,ai in a.items():
        for j,bj in b.items():
            if i+j<=n: out[i+j]=out.get(i+j,Fraction(0))+ai*bj
    return {k:v for k,v in out.items() if v}
def spow(a,p,n=ORD):
    out={0:Fraction(1)}
    for _ in range(p): out=smul(out,a,n)
    return out
sech={2*k:Fraction(sp.euler(2*k))/sp.factorial(2*k) for k in range(ORD//2+1)}
gd={2*k+1:sech[2*k]/(2*k+1) for k in range(ORD//2+1) if 2*k+1<=ORD}
uu={k:v for k,v in sech.items() if k>0}
sech12={0:Fraction(1)}
for j in range(1,ORD+1):
    c=Fraction(1)
    for i in range(j): c*=Fraction(1,2)-i
    c/=sp.factorial(j)
    for k,v in spow(uu,j).items(): sech12[k]=sech12.get(k,Fraction(0))+c*v
sech12={k:v for k,v in sech12.items() if v}
H32=smul({k:Fraction(1,2**k*sp.factorial(k)) for k in range(ORD+1)}, smul(sech,sech12))
e={}
for r in range(6):
    ss=smul(spow(gd,2*r),H32)
    for m in range(13): e[(r,m)]=ss.get(m,Fraction(0))
theirs={}
with open('../K3_fourier_all_orders/e_rm.pkl','rb') as f:
    for key,val in pickle.load(f).items():
        r,m=map(int,key.split(',')); theirs[(r,m)]=sp.Rational(val)
report("B1 e_{r,m} exact rebuild == K3 table (50 entries)", all(sp.Rational(e[k])==theirs[k] for k in theirs))
report("B2 e_{r,2r}=1, e_{r,2r+1}=1/2 (all r<=5)",
       all(e[(r,2*r)]==1 for r in range(6)) and all(e[(r,2*r+1)]==Fraction(1,2) for r in range(6)))
def eta(x): return (1-mp.power(2,1-x))*mp.zeta(x)
Ic={}
def I(m,j):
    if (m,j) not in Ic: Ic[(m,j)]=mp.quad(lambda ss: ss**m*mp.log1p(mp.e**(-ss))**j,[0,0.5,2,8,32,mp.inf])
    return Ic[(m,j)]
report("B3 I_{m,1} = m! eta(m+2) (m=0..5), I_{0,2} = zeta(3)/4",
       all(abs(I(m,1)-mp.factorial(m)*eta(m+2))<mp.mpf('1e-60') for m in range(6)) and abs(I(0,2)-mp.zeta(3)/4)<mp.mpf('1e-60'))
em={k: mp.mpf(str(v)) for k,v in e.items()}
def Dt(n,r): return 2*mp.fsum(em[(r,n-1-j)]*I(n-1-j,j)/mp.factorial(j) for j in range(1,n-2*r))
ev = all(abs(Dt(2*n,n-1)-2*mp.factorial(2*n-2)*eta(2*n))<mp.mpf('1e-55') for n in range(2,6))
od = all(abs(Dt(2*n+1,n-1)-mp.factorial(2*n-2)*(mp.zeta(2*n+1)+2*mp.fsum(eta(2*j)*mp.zeta(2*n+1-2*j) for j in range(1,n))))<mp.mpf('1e-55') for n in range(1,5))
report("B4 even tops D~_{2n,n-1} = 2(2n-2)! eta(2n) != 0, n=2..5", ev)
report("B5 odd tops D~_{2n+1,n-1} = (2n-2)![z(2n+1)+2 sum eta zeta] != 0, n=1..4", od)
lc = all(abs(Dt(2*m,m-1)*(-4)**(m-1)/mp.factorial(2*m-2) - (-1)**(m-1)*(2**(2*m-1)-1)*mp.zeta(2*m))<mp.mpf('1e-55') for m in range(1,7))
report("B6 leading coeff Theta_{2m} = (-1)^{m-1}(2^{2m-1}-1) zeta(2m), m=1..6", lc)
def Theta(n,u): return mp.fsum(Dt(n,m)*(-4*u)**m/mp.factorial(2*m) for m in range(0,n//2))
Q = pickle.load(open('../K3_fourier_all_orders/Qpolys.pkl','rb'))
u_sym=sp.Symbol('u'); zsub={sp.Symbol(f'z{k}'):sp.zeta(k) for k in range(2,13)}
okQ=True
for n in range(2,9):
    qn=sp.sympify(Q[n], locals={f'z{k}':sp.Symbol(f'z{k}') for k in range(2,13)})
    for uv in [0,1,2,5,7]:
        if abs(Theta(n,uv)-mp.mpf(str(qn.subs({u_sym:uv,**zsub}).evalf(50))))>mp.mpf('1e-40'): okQ=False
report("B7 Theta_n(u) [my D~-build] == K3 Qpolys (zeta products), n<=8, u in {0,1,2,5,7}", okQ)
# a_k quadrature extraction
def a_k(k,tau):
    tau=mp.mpf(tau); q=1/tau
    return mp.mpf(4)/mp.pi*mp.quad(lambda th:(mp.cos(th)**q+mp.sin(th)**q)**tau*mp.cos(4*k*th),[0,mp.pi/4,mp.pi/2])
taus=[mp.mpf(x) for x in ['0.001','0.002','0.003','0.004','0.006']]
Ak={k:-(-1)**k*4*mp.sqrt(2)/(mp.pi*(16*k*k-1)) for k in range(5)}
ak=[[a_k(k,t) for t in taus] for k in range(5)]
coef={}
for k in range(5):
    A=mp.matrix(len(taus),5); b=mp.matrix(len(taus),1)
    for i,tt in enumerate(taus):
        for j,n in enumerate(range(2,7)): A[i,j]=tt**n
        b[i]=ak[k][i]-Ak[k]
    sol=mp.lu_solve(A.T*A,A.T*b)
    coef[k]={n:(-1)**k*mp.pi/mp.sqrt(2)*sol[n-2] for n in range(2,7)}
thrA={2:'1e-8',3:'1e-5',4:'1e-4'}; okA=True
for n in range(2,5):
    mxs=max(abs(coef[k][n]-Theta(n,k*k))/max(abs(Theta(n,k*k)),mp.mpf('1e-30')) for k in range(5))
    print(f"   fresh a_k extraction: max rel.err n={n} over k=0..4 = {mp.nstr(mxs,2)} (threshold {thrA[n]})")
    okA &= mxs < mp.mpf(thrA[n])
report("B8a fresh a_k quadrature extraction matches Theta_n(k^2), n=2..4, k<=4", okA)
# full saved dataset (90 dps, 8 tau, k<=6): per-n thresholds matching asymptotic extraction error
dd=json.load(open('quad_ak.json')); tt=[mp.mpf(x) for x in dd['taus']]
AK6={k:-(-1)**k*4*mp.sqrt(2)/(mp.pi*(16*k*k-1)) for k in range(7)}
cf={}
for k in range(7):
    vs=[mp.mpf(x)-AK6[k] for x in dd['data'][str(k)]]
    A=mp.matrix(len(tt),7); b=mp.matrix(len(tt),1)
    for i,t0 in enumerate(tt):
        for j,n in enumerate(range(2,9)): A[i,j]=t0**n
        b[i]=vs[i]
    sol=mp.lu_solve(A.T*A,A.T*b)
    cf[k]={n:(-1)**k*mp.pi/mp.sqrt(2)*sol[n-2] for n in range(2,9)}
thr={2:'1e-9',3:'1e-6',4:'1e-5',5:'1e-3',6:'5e-3'}
okB=True
for n in range(2,7):
    mxs=max(abs(cf[k][n]-Theta(n,k*k))/max(abs(Theta(n,k*k)),mp.mpf('1e-30')) for k in range(7))
    print(f"   saved dataset: max rel.err n={n} over k=0..6 = {mp.nstr(mxs,2)} (threshold {thr[n]})")
    okB &= mxs < mp.mpf(thr[n])
report("B8b saved 90-dps dataset: extraction matches Theta_n, n<=6, k<=6", okB)
# degree overfit test on extracted values
degok=True
for n in range(2,7):
    deg=n//2-1; us=[k*k for k in range(5)]; vs=[coef[k][n] for k in range(5)]
    A=mp.matrix(5,deg+2); b=mp.matrix(5,1)
    for i,u in enumerate(us):
        for j in range(deg+2): A[i,j]=mp.mpf(u)**j
        b[i]=vs[i]
    sol=mp.lu_solve(A.T*A,A.T*b)
    scale=max(abs(v) for v in vs)
    if abs(sol[deg+1]*16**(deg+1))/scale > mp.mpf('0.02'): degok=False
report("B9 overfit u^{deg+1} coefficient negligible (degree exactly floor(n/2)-1), n<=6", degok)

print("="*80)
fails=[n for n,ok in PASS if not ok]
print(f"SUMMARY: {len(PASS)-len(fails)}/{len(PASS)} checks PASS; failures: {fails if fails else 'NONE'}; runtime {time.time()-t00:.0f}s")
