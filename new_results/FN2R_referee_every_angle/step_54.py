verify = r'''#!/usr/bin/env python
"""verify_fn2r.py — consolidated FN2-referee verification (independent of FN2's verify_fn2.py).
  R1  Lemma B: d_l(c)=(-1)^l(1-c^2)P'_{l-1}(c)/l, l<=10; c=0 -> -C(-1/2,m)   [sympy exact]
  R2  theta_{ij} (own reversion) + master formula -> displayed J1,J2,J3, general c [sympy exact]
  R3  paper thm:P2 exact series (f=x) through eps^3 == displayed J's, general c    [sympy exact]
  R4  cumulants k1,k2,k3(theta) vs closed forms; Mellin=Legendre at w=+-0.9        [numeric]
  R5  residual after J3 scaled by eps^4 log^4 n bounded (exp, linexp; pi/3, 2pi/3) [numeric]
  R6  resid = eps^4 J4 (1+o(1)) at one point                                       [numeric]
  R7  contact sin(pi x/2): D/L0, D/L3, D/full-model, refined-law scaled residual   [numeric]
  R8  gap regime m=3/2 and m=1/2 exact-Phi models: |1-D/Fm|                        [numeric]
  R9  Klayer = (1-c) k1 / 2                                                        [numeric]
Target < 5 min. All numbers printed.
"""
import time
t_start = time.time()
def el(): return "%5.1fs" % (time.time()-t_start)

# ---------------- R1-R3 symbolic ----------------
import sympy as sp
c, z, a, A, B = sp.symbols('c z a A B')
t = sp.Symbol('t')
print("== R1: Lemma B tail coefficients ==")
phi_t = 1 - (1+c*t)*(1+2*c*t+t**2)**sp.Rational(-1,2)
ser = sp.series(phi_t, t, 0, 11).removeO().expand()
ok1 = all(sp.simplify(ser.coeff(t,l) - (-1)**l*(1-c**2)*sp.legendre(l-1,c).diff(c)/l)==0 for l in range(1,11))
ok1b = all(sp.simplify(((-1)**(2*m)*(1-c**2)*sp.legendre(2*m-1,c).diff(c)/(2*m)).subs(c,0)
           + sp.binomial(sp.Rational(-1,2),m))==0 for m in range(1,5))
ok1c = all(sp.simplify(k*(sp.legendre(k,sp.Symbol('x')) - sp.Symbol('x')*sp.legendre(k-1,sp.Symbol('x')))
           - (sp.Symbol('x')**2-1)*sp.legendre(k-1,sp.Symbol('x')).diff(sp.Symbol('x')))==0 for k in range(1,7))
print("  d_l l<=10:", ok1, " c=0 -> -C(-1/2,m):", ok1b, " Legendre identity:", ok1c, " [%s]"%el())

print("== R2: theta_ij reversion + master -> displayed ==")
t, W, mu = sp.symbols('t W mu')
K = 3
u = -a*t + A*t**2/2 - B*t**3/6
Acal = sp.expand(sum((-1)**r*u**r/r for r in range(1,K+2)))
Acal = sum(mm for mm in Acal.args if sp.degree(mm,t) <= K+1)
v = (-A*t + B*t**2/2)/a
Btil = sp.expand(sum((-1)**(r+1)*v**r/r for r in range(1,K+1)))
Btil = sum(mm for mm in Btil.args if sp.degree(mm,t) <= K)
idx = [(i,j) for i in range(1,K+1) for j in range(0,K+1-i)]
vs = sp.symbols(' '.join('t%d%d'%(i,j) for i,j in idx))
Tau = sum(vv*W**i*mu**j for vv,(i,j) in zip(vs,idx))
def trunc(e):
    e = sp.expand(e); return sum(mm for mm in e.args if sp.degree(mm,W)+sp.degree(mm,mu) <= K)
eq = trunc(trunc(Acal.subs(t,Tau)) - mu*trunc(Btil.subs(t,Tau)) - W)
sols = {}
for (i,j),vv in zip(idx,vs):
    sol = sp.solve(sp.expand(eq.subs(sols)).coeff(W,i).coeff(mu,j), vv)
    sols[vv] = sol[0]
th = {k2_: sp.simplify(sols[vv]) for k2_,vv in zip(idx,vs)}
k2, k3, Xi = sp.symbols('k2 k3 Xi')
mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
M = {1: 2*z*mu1, 2: 2*z*mu2, 3: 2*z*mu3}
d = lambda l: (-1)**l*(1-c**2)*sp.legendre(l-1, c).diff(c)/l
Tser = 0
for (i,j),tij in th.items():
    for l in range(2, K+2):
        Tser += sp.series(tij*mu**(i+j)*(-1)**(i+1)*sp.factorial(i)*d(l)/(l-1)**(i+1)*(mu/(a*(1+mu)))**(l-1), mu, 0, K+1).removeO()
T2 = sp.expand(Tser).coeff(mu,2); T3 = sp.expand(Tser).coeff(mu,3)
J1 = sp.expand(th[(1,0)]*M[1])
J2 = sp.expand(th[(2,0)]*M[2] + th[(1,1)]*M[1] + T2)
J3 = sp.expand(th[(3,0)]*M[3] + th[(2,1)]*M[2] + th[(1,2)]*M[1] + T3)
J1p = 2*z*mu1/a
J2p = z*(A*(mu2-2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
J3p = ( z*(a**4*mu3 - (3*A*a**2+B*a)*(mu3-3*mu2) + 3*A**2*(mu3-4*mu2+2*mu1))/(3*a**5)
        - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
sub = {z: (1-c)/2}
ok2 = (sp.simplify((J1-J1p).subs(sub))==0 and sp.simplify(sp.expand((J2-J2p).subs(sub)))==0
       and sp.simplify(sp.expand((J3-J3p).subs(sub)))==0)
Q2 = -(a**2-A)/(2*a**3)*(Xi**2+k2) - A*Xi/a**3 + 1/(2*a**2)
lam3 = a**4-3*a**2*A-a*B+3*A**2
Q3 = lam3/(6*a**5)*mu3 + (3*a**2*A+a*B-4*A**2)/(2*a**5)*mu2 + A**2*Xi/a**5 - A/a**4
subR = {z: sp.Rational(1,2), c: 0}
ok2b = (sp.simplify((J1-Xi/a).subs(subR))==0 and sp.simplify(sp.expand((J2-Q2).subs(subR)))==0
        and sp.simplify(sp.expand((J3-Q3).subs(subR)))==0)
print("  master - displayed at general c = 0:", ok2, " ; pi/2 -> P3 Q's = 0:", ok2b, " [%s]"%el())

print("== R3: thm:P2 exact series (f=x) through eps^3 ==")
eps, L = sp.symbols('eps L')
a1, a2, a3 = sp.symbols('a1 a2 a3')
w = sp.Symbol('w')
logLam = a1*w + (a2 - sp.pi**2/6)*w**2 + a3*w**3
Lam = sp.expand(sp.exp(logLam).series(w,0,4).removeO())
al1 = Lam.coeff(w,1); al2 = Lam.coeff(w,2); al3 = Lam.coeff(w,3)
n_min_eps = sp.expand(sp.exp(-eps*L).series(eps,0,4).removeO())
mainfac = sp.expand((1 + (sp.pi**2/6)*eps**2) * (1 + al1*eps + al2*eps**2 + al3*eps**3) * n_min_eps)
term1 = sp.expand(2*z - 2*z*mainfac)
term1 = sum(mm for mm in term1.args if sp.degree(mm,eps)<=3)
Lam1 = (sp.legendre(0,c)-sp.legendre(2,c))/3
Lam2 = (sp.legendre(1,c)-sp.legendre(3,c))/5
mpart = ( sp.series(eps/(1-eps), eps, 0, 4).removeO()*Lam1*sp.series(eps/(1+eps), eps, 0, 4).removeO()
        - sp.series(eps/(2-eps), eps, 0, 3).removeO()*Lam2*sp.series(eps**2/(1+eps)**2, eps, 0, 3).removeO() )
two_minus_L = sp.expand(sum(mm for mm in sp.expand(term1+mpart).args if sp.degree(mm,eps)<=3))
k1s, k2s, k3s = -a1, 2*a2, -6*a3
Xis = L + k1s
J1x = 2*z*Xis; J2x = z*(2*(1-z) - (Xis**2+k2s)); J3x = z*(Xis**3+3*k2s*Xis+k3s)/3 - c*z*(1-z)
ok3 = all(sp.simplify(sp.expand((two_minus_L.coeff(eps,i) - [J1x,J2x,J3x][i-1]).subs(sub)))==0 for i in (1,2,3))
print("  exact series - displayed J (f=x), general c = 0:", ok3, " [%s]"%el())

# ---------------- numeric ----------------
import mpmath as mp
mp.mp.dps = 50
def hc(s, c_):
    if s <= 0: return mp.mpf(0)
    r = mp.sqrt(1 + 2*c_*s + s**2)
    return (1-c_) - (1-c_**2)/(r + s + c_)
def phi_s(s,c_):
    r = mp.sqrt(1+2*c_*s+s**2)
    return (1-c_**2)/(r*(r+s+c_))
def S12(zv):
    return mp.zeta(3) - mp.polylog(3, 1-zv) + mp.log(1-zv)*mp.polylog(2, 1-zv) \
           + mp.mpf('0.5')*mp.log(1-zv)**2*mp.log(zv)
def kforms(thv):
    zv = mp.sin(thv/2)**2; w2 = 1/mp.tan(thv/2)**2; u = w2*mp.log(mp.cos(thv/2)**2)
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,zv)-u**2, 2+2*u**3-12*w2*S12(zv))

print("== R4: cumulants + Mellin ==")
worst = mp.mpf(0)
for thv in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    cv = mp.cos(thv)
    raw=[]
    for r in (1,2,3):
        i1 = mp.quad(lambda u: phi_s(mp.exp(-u),cv)*mp.exp(-u)*u**r, [0,mp.inf])
        i2 = mp.quad(lambda u: phi_s(mp.exp(u),cv)*mp.exp(u)*(-u)**r, [0,mp.inf])
        raw.append((i1+i2)/(1-cv))
    m1,m2,m3 = raw
    kq = (m1, m2-m1**2, m3-3*m1*m2+2*m1**3)
    kf = kforms(thv)
    worst = max(worst, max(abs(kq[i]-kf[i]) for i in range(3)))
cv = mp.cos(mp.pi/3); zv = mp.sin(mp.pi/6)**2
p = 10
worstM = mp.mpf(0)
for wv in [mp.mpf('-0.9'), mp.mpf('0.9')]:
    i1 = mp.quad(lambda u: (u**p)**wv*phi_s(u**p,cv)*p*u**(p-1), [0, 1])
    i2 = mp.quad(lambda u: (u**-p)**wv*phi_s(u**-p,cv)*p*u**(-p-1), [0, 1])
    lhs = (i1+i2)/(1-cv)
    rhs = (mp.pi*wv/mp.sin(mp.pi*wv))*2*zv*mp.hyper([-wv,1+wv],[2],zv)/(1-cv)
    worstM = max(worstM, abs(lhs-rhs))
print("  cumulants worst |diff| (3 angles):", mp.nstr(worst,2), " Mellin=Legendre w=+-0.9:", mp.nstr(worstM,2), " [%s]"%el())

print("== R5: residual after J3 / eps^4 log^4 n ==")
def Jpaper(thv, a_, A_, B_, n):
    k1,k2,k3 = kforms(thv)
    zv = mp.sin(thv/2)**2; cv = mp.cos(thv)
    Xi = mp.log(n*a_) + k1
    mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
    return (2*zv*mu1/a_,
            zv*(A_*(mu2-2*mu1) - a_**2*mu2 + 2*a_*(1-zv))/a_**3,
            ( zv*(a_**4*mu3 - (3*A_*a_**2+B_*a_)*(mu3-3*mu2) + 3*A_**2*(mu3-4*mu2+2*mu1))/(3*a_**5)
              - 4*A_*zv*(1-zv)/a_**4 - cv*zv*(1-zv)/a_**3 ))
def make(which):
    if which=='exp':
        return (lambda x: (mp.exp(x)-1)/(mp.e-1), lambda x: mp.exp(x)/(mp.e-1),
                (mp.e/(mp.e-1),)*3)
    return (lambda x: (1+x)/2*mp.exp(x-1), lambda x: (2+x)/2*mp.exp(x-1),
            (mp.mpf('1.5'), mp.mpf(2), mp.mpf('2.5')))
def Dn(which, thv, n):
    f, fp, _ = make(which)
    nn = mp.mpf(n); cv = mp.cos(thv); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        fx = f(1-tt)
        if fx <= 0: return mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fp(1-tt))
        return hc(mp.exp(lyp), cv)*mp.exp(v)/nn
    pts = sorted(set([-80, -20, lnn-8, lnn-2, lnn+2, lnn]))
    return mp.quad(bump, pts, maxdegree=9)
mx = mp.mpf(0)
for which in ['exp','linexp']:
    a_,A_,B_ = make(which)[2]
    for thn, thv in [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3)]:
        for n in [10**4, 10**5, 10**6]:
            D = Dn(which, thv, n)
            J1,J2,J3 = Jpaper(thv, a_, A_, B_, mp.mpf(n))
            eps_ = 1/mp.mpf(n-1)
            sc = (D - (J1*eps_+J2*eps_**2+J3*eps_**3))/(eps_**4*mp.log(n)**4)
            mx = max(mx, abs(sc))
            print("  f=%6s th=%5s n=%7d: scaled resid = %s" % (which, thn, n, mp.nstr(sc,6)))
print("  max |scaled| =", mp.nstr(mx,4), " [%s]"%el())

print("== R6: resid = eps^4 J4 (1+o(1)) at exp, pi/3, n=1e5 ==")
thv = mp.pi/3; cv = mp.cos(thv)
raw=[]
for r in (1,2,3,4):
    i1 = mp.quad(lambda u: phi_s(mp.exp(-u),cv)*mp.exp(-u)*u**r, [0,mp.inf])
    i2 = mp.quad(lambda u: phi_s(mp.exp(u),cv)*mp.exp(u)*(-u)**r, [0,mp.inf])
    raw.append((i1+i2)/(1-cv))
m1,m2,m3,m4 = raw
k4q = m4 - 4*m3*m1 - 3*m2**2 + 12*m2*m1**2 - 6*m1**4
a_,A_,B_ = (mp.e/(mp.e-1),)*3; D4 = a_
n = 10**5
k1,k2,k3 = kforms(thv); zv = mp.sin(thv/2)**2
Xi = mp.log(n*a_) + k1
mu4 = Xi**4 + 6*k2*Xi**2 + 4*k3*Xi + (k4q+3*k2**2)
M1,M2,M3,M4 = 2*zv*Xi, 2*zv*(Xi**2+k2), 2*zv*(Xi**3+3*k2*Xi+k3), 2*zv*mu4
th40 = (15*A_**3 - 10*A_*B_*a_ + 7*A_*a_**4 + 6*B_*a_**3 - a_**6 + a_**2*(-18*A_**2 + D4))/(24*a_**7)
th31 = -(7*a_**4*A_+6*a_**3*B_-24*a_**2*A_**2+a_**2*D4-13*a_*A_*B_+23*A_**3)/(6*a_**7)
th22 = -3*A_*(2*a_**2*A_+a_*B_-3*A_**2)/(2*a_**7)
th13 = -A_**3/a_**7
T4 = -(cv-1)*(cv+1)*(96*A_**2 - 48*A_*a_**2 + 9*A_*a_*cv - 24*B_*a_ + 12*a_**4 + 9*a_**3*cv + 5*a_**2*cv**2 - a_**2)/(24*a_**6)
J4 = th40*M4 + th31*M3 + th22*M2 + th13*M1 + T4
D = Dn('exp', thv, n)
J1,J2,J3 = Jpaper(thv, a_, A_, B_, mp.mpf(n))
eps_ = 1/mp.mpf(n-1)
resid = D - (J1*eps_+J2*eps_**2+J3*eps_**3)
print("  resid/(eps^4 J4) =", mp.nstr(resid/(eps_**4*J4),8), " [%s]"%el())

print("== R7: contact sin(pi x/2) ==")
beta = mp.pi**2/8; m = mp.mpf(2); cm = 1-1/m; ym = cm*(1-mp.log(cm))
def Uroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, -1)
def Vroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, 0)
def ell(y):
    if y <= ym: return mp.mpf(0)
    return (Uroot(y)**(1/m) - Vroot(y)**(1/m)).real
def Ffun(y): return Uroot(y)**(1/m)
def F2f(y):
    U = Uroot(y); Up = U/(U-cm); Upp = -cm*Up/(U-cm)**2
    return (1/m)*((1/m-1)*U**(1/m-2)*Up**2 + U**(1/m-1)*Upp)
fs  = lambda x: mp.sin(mp.pi*x/2)
fps = lambda x: mp.pi/2*mp.cos(mp.pi*x/2)
def Dn_sin(thv, n):
    nn = mp.mpf(n); cv = mp.cos(thv); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        fx = fs(1-tt)
        if fx <= 0: return mp.mpf(0)
        lyp = mp.log(nn) + (nn-1)*mp.log(fx) + mp.log(fps(1-tt))
        return hc(mp.exp(lyp), cv)*mp.exp(v)/nn
    mid = (lnn + mp.log(lnn))/2
    return mp.quad(bump, sorted(set([-60, mid-8, mid-2, mid+2, mid+8, lnn])), maxdegree=9)
def Fm(thv, n):
    cv = mp.cos(thv); Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn); us = Lam - ym
    I = mp.quad(lambda u: ell(Lam-u)*phi_s(mp.exp(u),cv)*mp.exp(u),
                sorted(set([-80, 0, us-4, us, us+2, us+8, 60])), maxdegree=9)
    return I/Nn
for thn, thv in [('pi/3', mp.pi/3), ('2pi/3', 2*mp.pi/3)]:
    cv = mp.cos(thv); k1,k2,k3 = kforms(thv)
    for n in [10**5, 10**6]:
        D = Dn_sin(thv, n)
        Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
        L0 = (1-cv)*(mp.log(n)/(m*beta*n))**(1/m)
        L3 = (1-cv)*(Ffun(Lam+k1) + (k2/2)*F2f(Lam+k1))/Nn
        Fmv = Fm(thv, n)
        pred = 1 + ((m-1)*mp.log(Lam) + m*k1)/(m**2*Lam)
        sc = (D/L0 - pred)*(Lam/mp.log(Lam))**2
        print("  th=%5s n=%8d: D/L0=%s D/L3=%s D/Fm=%s scaled(D/L0-pred)=%s" % (
            thn, n, mp.nstr(D/L0,7), mp.nstr(D/L3,7), mp.nstr(D/Fmv,9), mp.nstr(sc,5)))
print("  [%s]"%el())

print("== R8: gap regime m=3/2 and m=1/2 (exact-Phi) ==")
def setup2(mv):
    m_ = mp.mpf(mv); cm_ = 1-1/m_
    if m_ > 1:
        ym_ = cm_*(1-mp.log(cm_))
        U = lambda y: -cm_*mp.lambertw(-mp.exp(-y/cm_)/cm_, -1)
        V = lambda y: -cm_*mp.lambertw(-mp.exp(-y/cm_)/cm_, 0)
        def ell_(y):
            if y <= ym_: return mp.mpf(0)
            return (U(y)**(1/m_) - V(y)**(1/m_)).real
    else:
        ym_ = -mp.inf
        U = lambda y: abs(cm_)*mp.lambertw(mp.exp(y/abs(cm_))/abs(cm_), 0)
        ell_ = lambda y: (U(y)**(1/m_)).real
    return m_, cm_, ym_, ell_
def Dn_exact(thv, n, m_, beta_):
    nn = mp.mpf(n); cv = mp.cos(thv); lnn = mp.log(nn)
    def bump(v):
        tt = mp.exp(v)/nn
        if tt > 1: return mp.mpf(0)
        lyp = mp.log(nn) + mp.log(beta_*m_*tt**(m_-1)) - nn*beta_*tt**m_
        return hc(mp.exp(lyp), cv)*mp.exp(v)/nn
    mid = (1-1/m_)*lnn + mp.log(lnn/beta_)/m_
    return mp.quad(bump, sorted(set([-80, mid-10, mid-3, mid+3, mid+10, lnn])), maxdegree=9)
thv = mp.pi/3; cv = mp.cos(thv)
for mv in [mp.mpf('1.5'), mp.mpf('0.5')]:
    m_, cm_, ym_, ell_ = setup2(mv)
    n = 10**5
    D = Dn_exact(thv, n, m_, mp.mpf(1))
    Nn = (mp.mpf(n))**(1/m_); Lam = mp.log(m_*Nn)
    if ym_ == -mp.inf: pts = [-80, 0, Lam-4, Lam, Lam+4, 60]
    else:
        us = Lam - ym_; pts = sorted(set([-80, 0, us-4, us, us+2, us+8, 60]))
    I = mp.quad(lambda u: ell_(Lam-u)*phi_s(mp.exp(u),cv)*mp.exp(u), pts, maxdegree=9)
    print("  m=%s n=1e5: |1 - D/Fm| = %s" % (mp.nstr(m_,4), mp.nstr(abs(1-D/(I/Nn)),3)))
print("  [%s]"%el())

print("== R9: Klayer identity ==")
thv = mp.pi/3; cv = mp.cos(thv); k1 = kforms(thv)[0]
i1 = mp.quad(lambda s: hc(s,cv)/s, [0, 1])
i2 = mp.quad(lambda s: (mp.sqrt(1+2*cv*s+s**2)-s-cv)/s, [1, mp.inf])
Kl = (i1-i2)/2
print("  Klayer - (1-c)k1/2 =", mp.nstr(Kl - (1-cv)*k1/2, 2), " [%s]"%el())

print("TOTAL runtime: %.1fs" % (time.time()-t_start))
'''
open("verify_fn2r.py","w").write(verify)
print("verify_fn2r.py written")
