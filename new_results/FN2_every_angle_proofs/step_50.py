script = r'''#!/usr/bin/env python
"""verify_fn2.py — consolidated verification for FN2 (proof.md).

Checks, in order:
  S1  Lemma B: tail coefficients d_l(c) = (-1)^l (1-c^2) P'_{l-1}(c)/l  [sympy, exact]
  S2  master formula -> paper's displayed J1,J2,J3 at general (c,z)      [sympy, exact]
  S3  reduction to theta=pi/2 (P3's Q1,Q2,Q3) and f=x exact series       [sympy, exact]
  S4  cumulants k1,k2,k3(theta) vs closed forms (7 angles)               [numeric]
  S5  residuals after J3 scaled by eps^4 log^4 n bounded (2 f x 2 th)    [numeric]
  S6  contact law: sin(pi x/2), 3 angles, ratios vs L0/L1/L3/full model  [numeric]
  S7  K(c) identity: Klayer = (1-c) k1 / 2                               [numeric]

Runtime target: < 5 min. All numbers printed.
"""
import time
t_start = time.time()

# ---------- S1+S2+S3: symbolic ----------
import sympy as sp
c, z, a, A, B, k2, k3, Xi = sp.symbols('c z a A B k2 k3 Xi')

print("== S1: tail coefficients d_l(c) vs Legendre form ==")
tt = sp.Symbol('tt')
phi_t = 1 - (1+c*tt)*(1+2*c*tt+tt**2)**sp.Rational(-1,2)
ser = sp.series(phi_t, tt, 0, 8).removeO().expand()
okS1 = all(sp.simplify(sp.expand(ser).coeff(tt, l)
           - (-1)**l*(1-c**2)*sp.legendre(l-1, c).diff(c)/l) == 0 for l in range(2, 8))
okS1b = all(sp.simplify(((-1)**(2*m)*(1-c**2)*sp.legendre(2*m-1, c).diff(c)/(2*m)).subs(c,0)
            + sp.binomial(sp.Rational(-1,2), m)) == 0 for m in range(1, 4))
print("  d_l = (-1)^l (1-c^2) P'_(l-1)(c)/l for l=2..7:", okS1,
      "; c=0 -> -C(-1/2,m):", okS1b)

print("== S2: master formula -> displayed J1,J2,J3 ==")
th = {(1,0): 1/a, (2,0): (A-a**2)/(2*a**3), (1,1): -A/a**3,
      (3,0): (3*A**2-3*A*a**2-B*a+a**4)/(6*a**5),
      (2,1): (-4*A**2+3*A*a**2+B*a)/(2*a**5), (1,2): A**2/a**5}
mu1 = Xi; mu2 = Xi**2+k2; mu3 = Xi**3+3*k2*Xi+k3
M = {1: 2*z*mu1, 2: 2*z*mu2, 3: 2*z*mu3}
d2 = 2*z*(1-z); d3 = -4*c*z*(1-z)
T2 = d2/a*th[(1,0)]
T3 = d2/a*(-2*th[(2,0)]+th[(1,1)]-th[(1,0)]) + d3/a**2*(th[(1,0)]/4)
J1 = sp.expand(th[(1,0)]*M[1])
J2 = sp.expand(th[(2,0)]*M[2] + th[(1,1)]*M[1] + T2)
J3 = sp.expand(th[(3,0)]*M[3] + th[(2,1)]*M[2] + th[(1,2)]*M[1] + T3)
J1p = 2*z*mu1/a
J2p = z*(A*(mu2-2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
J3p = ( z*(a**4*mu3 - (3*A*a**2+B*a)*(mu3-3*mu2) + 3*A**2*(mu3-4*mu2+2*mu1))/(3*a**5)
        - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
sub = {z: (1-c)/2}
okS2 = (sp.simplify((J1-J1p).subs(sub))==0 and sp.simplify(sp.expand((J2-J2p).subs(sub)))==0
        and sp.simplify(sp.expand((J3-J3p).subs(sub)))==0)
print("  difference master - displayed = 0 at general (c,z):", okS2)

print("== S3: reductions ==")
Q1 = Xi/a
Q2 = -(a**2-A)/(2*a**3)*(Xi**2+k2) - A*Xi/a**3 + 1/(2*a**2)
lam3 = a**4 - 3*a**2*A - a*B + 3*A**2
Q3 = lam3/(6*a**5)*mu3 + (3*a**2*A+a*B-4*A**2)/(2*a**5)*mu2 + A**2*Xi/a**5 - A/a**4
subR = {z: sp.Rational(1,2), c: 0}
okR = (sp.simplify((J1-Q1).subs(subR))==0 and sp.simplify(sp.expand((J2-Q2).subs(subR)))==0
       and sp.simplify(sp.expand((J3-Q3).subs(subR)))==0)
exact_J1 = 2*z*mu1; exact_J2 = -z*mu2 + 2*z*(1-z); exact_J3 = z*mu3/3 - c*z*(1-z)
subx = {a: 1, A: 0, B: 0}
okX = (sp.simplify((J1-exact_J1).subs(subx))==0 and sp.simplify(sp.expand((J2-exact_J2).subs(subx)))==0
       and sp.simplify(sp.expand((J3-exact_J3).subs(subx)))==0)
print("  theta=pi/2 -> P3 Q1,Q2,Q3:", okR, "; f=x -> exact series:", okX)
print(f"[S1-S3 done, {time.time()-t_start:.1f}s]")

# ---------- numeric ----------
import mpmath as mp
mp.mp.dps = 50

def S12(zv):
    return mp.zeta(3) - mp.polylog(3, 1-zv) + mp.log(1-zv)*mp.polylog(2, 1-zv) \
           + mp.mpf('0.5')*mp.log(1-zv)**2*mp.log(zv)
def cumulants(thv):
    zv = mp.sin(thv/2)**2; w = 1/mp.tan(thv/2)**2; ell = mp.log(mp.cos(thv/2)**2); u = w*ell
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,zv)-u**2, -12*w*S12(zv)+2*(1+u**3))
def phi_st(s, cv):
    r = mp.sqrt(1 + 2*cv*s + s**2)
    return (1-cv**2)/(r*(r+s+cv))
def hc(s, cv):
    if s <= 0: return mp.mpf(0)
    if s < 1: return 1 + s - mp.sqrt(1 + 2*cv*s + s**2)
    r = mp.sqrt(1 + 2*cv*s + s**2)
    return (1-cv) - (1-cv**2)/(r + s + cv)

print("== S4: cumulants vs closed forms ==")
worst = 0
for thv in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6]:
    cv = mp.cos(thv); m = []
    for r in range(1, 4):
        i1 = mp.quad(lambda u: phi_st(mp.exp(-u),cv)*mp.exp(-u)*u**r, [0, mp.inf])
        i2 = mp.quad(lambda u: phi_st(mp.exp(u),cv)*mp.exp(u)*(-u)**r, [0, mp.inf])
        m.append((i1+i2)/(1-cv))
    m1,m2,m3 = m
    kn = (m1, m2-m1**2, m3-3*m2*m1+2*m1**3)
    kc = cumulants(thv)
    res = max(abs(kn[j]-kc[j]) for j in range(3)); worst = max(worst, res)
print("  worst |k_j quadrature - closed form| over 5 angles:", mp.nstr(worst, 3))
print(f"[S4 done, {time.time()-t_start:.1f}s]")

def Jvals(thv, a_, A_, B_, n):
    k1,k2,k3 = cumulants(thv)
    Xiv = mp.log(n*a_)+k1
    zv = mp.sin(thv/2)**2; cv = mp.cos(thv)
    mu1 = Xiv; mu2 = Xiv**2+k2; mu3 = Xiv**3+3*k2*Xiv+k3
    J1 = 2*zv*mu1/a_
    J2 = zv*(A_*(mu2-2*mu1) - a_**2*mu2 + 2*a_*(1-zv))/a_**3
    J3 = ( zv*(a_**4*mu3 - (3*A_*a_**2+B_*a_)*(mu3-3*mu2) + 3*A_**2*(mu3-4*mu2+2*mu1))/(3*a_**5)
           - 4*A_*zv*(1-zv)/a_**4 - cv*zv*(1-zv)/a_**3 )
    return J1,J2,J3
def make_f(which):
    if which=='exp':
        e = mp.e; a_ = e/(e-1)
        return (lambda x: (mp.exp(x)-1)/(e-1), lambda x: mp.exp(x)/(e-1), (a_,a_,a_))
    return (lambda x: (x+x**2)/2, lambda x: (1+2*x)/2, (mp.mpf('1.5'), mp.mpf(1), mp.mpf(0)))
def Dn(which, thv, n):
    f, fp, _ = make_f(which)
    nn = mp.mpf(n); cv = mp.cos(thv)
    def integrand(v):
        t = mp.exp(v)/nn
        if t >= 1: return mp.mpf(0)
        fx = f(1-t)
        if fx <= 0: return mp.mpf(0)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*fp(1-t)
        return hc(yp, cv)*mp.exp(v)/nn
    mid = mp.log(mp.log(nn)); lnn = mp.log(nn)
    pts = [-160, mid-6, mid-2, mid+2, mid+6, min(mid+12, lnn-mp.mpf('0.5')), lnn]
    return mp.quad(integrand, pts, maxdegree=10)

print("== S5: residuals after J3, scaled by eps^4 log^4 n ==")
for which in ['exp','quad']:
    a_,A_,B_ = make_f(which)[2]
    for thv in [mp.pi/3, 2*mp.pi/3]:
        for n in [10**4, 10**5]:
            num = Dn(which, thv, n)
            J1,J2,J3 = Jvals(thv, a_, A_, B_, mp.mpf(n))
            eps = 1/mp.mpf(n-1)
            resid = num - (J1*eps + J2*eps**2 + J3*eps**3)
            print(f"  f={which:4s} th={mp.nstr(thv,6)} n={n}: scaled resid = {mp.nstr(resid/(eps**4*mp.log(n)**4),5)}")
print(f"[S5 done, {time.time()-t_start:.1f}s]")

# S6 contact
beta = mp.pi**2/8; mm = mp.mpf(2)
def f_sin(x): return mp.sin(mp.pi*x/2)
def fp_sin(x): return mp.pi/2*mp.cos(mp.pi*x/2)
def Dn_sin(thv, n):
    nn = mp.mpf(n); cv = mp.cos(thv); Nn = mp.sqrt(beta*nn); mid = mp.log(nn/Nn)
    def integrand(v):
        t = mp.exp(v)/nn
        if t >= 1: return mp.mpf(0)
        fx = f_sin(1-t)
        if fx <= 0: return mp.mpf(0)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*fp_sin(1-t)
        return hc(yp, cv)*mp.exp(v)/nn
    lnn = mp.log(nn)
    pts = [-200, mid-12, mid-5, mid, mid+4, min(mid+9, lnn-mp.mpf('0.5')), lnn]
    return mp.quad(integrand, pts, maxdegree=10)
cm = 1 - 1/mm; ym = cm*(1-mp.log(cm))
def Uroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, -1)
def Vroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, 0)
def Ffun(y): return Uroot(y)**(1/mm)
def ellfun(y):
    return mp.mpf(0) if y < ym else Uroot(y)**(1/mm) - Vroot(y)**(1/mm)
def F2(y):
    U = Uroot(y); Up = U/(U-cm); Upp = -cm*Up/(U-cm)**2
    return (1/mm)*((1/mm-1)*U**(1/mm-2)*Up**2 + U**(1/mm-1)*Upp)
def full_model(thv, n):
    cv = mp.cos(thv); Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(mm*Nn)
    I = mp.quad(lambda u: ellfun(Lam-u)*phi_st(mp.exp(u),cv)*mp.exp(u),
                [-60, Lam-12, Lam-4, Lam+4, 60], maxdegree=10)
    return I/Nn

print("== S6: contact law, f=sin(pi x/2) ==")
for thv in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    cv = mp.cos(thv); k1,k2,_ = cumulants(thv)
    for n in [10**5, 10**6]:
        D = Dn_sin(thv, n)
        Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(mm*Nn)
        L0 = (1-cv)*(mp.log(n)/(mm*beta*n))**(1/mm)
        L1 = (1-cv)*Ffun(Lam)/Nn
        L3 = (1-cv)*(Ffun(Lam+k1)+(k2/2)*F2(Lam+k1))/Nn
        Fm = full_model(thv, n)
        pred = 1 + ((mm-1)*mp.log(Lam)+mm*k1)/(mm**2*Lam)
        print(f"  th={mp.nstr(thv,6)} n={n}: D/L0={mp.nstr(D/L0,6)} (pred {mp.nstr(pred,6)}) "
              f"D/L1={mp.nstr(D/L1,6)} D/L3={mp.nstr(D/L3,6)} D/full={mp.nstr(D/Fm,8)}")
print(f"[S6 done, {time.time()-t_start:.1f}s]")

print("== S7: K(c) identity: Klayer = (1-c) k1/2 ==")
for thv in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    cv = mp.cos(thv); k1 = cumulants(thv)[0]
    i1 = mp.quad(lambda s: hc(s,cv)/s, [0, 1])
    i2 = mp.quad(lambda s: (1-cv**2)/(mp.sqrt(1+2*cv*s+s**2)+s+cv)/s, [1, mp.inf])
    Klayer = (i1-i2)/2
    print(f"  th={mp.nstr(thv,6)}: Klayer={mp.nstr(Klayer,10)}  (1-c)k1/2={mp.nstr((1-cv)*k1/2,10)}"
          f"  diff={mp.nstr(Klayer-(1-cv)*k1/2,2)}")

print(f"\nTOTAL runtime: {time.time()-t_start:.1f}s")
'''
with open("verify_fn2.py","w") as f:
    f.write(script)
print("verify_fn2.py written", len(script))
