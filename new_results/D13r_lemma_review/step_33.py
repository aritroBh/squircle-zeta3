import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
yJ = mp.sqrt(pi/(2*mp.mpf(10)**10))
mn = mp.mpf('9'); mn2 = mp.mpf('9'); mxB = mp.mpf('-9')
random.seed(46)
for _ in range(3000):
    y = yJ*mp.mpf(10)**random.uniform(-6, 0) if random.random()<0.5 else yJ*random.random()
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    sig = (rho-1)/y**2; x0=1/(rho*y); t=rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t)
    Ap = 1 - mp.atan(zz)/zz
    breg_h = -sig/(2*rho) + (-2*y*sig + 2*sqrt2*sig + sqrt2*rho - 2*sqrt2*rho*Ap)/(2*sqrt2*(2+sqrt2*t))
    Btrue = (phi_true + pi/4)/((pi/sqrt2)*y)
    mn = min(mn, breg_h - Btrue)          # need >= 0 (given both positive)
    mn2 = min(mn2, Btrue)                 # need > 0 for |Btrue|=Btrue
    mxB = max(mxB, Btrue)
print("final cell true pairs: min(breg_h[A1'] - Btrue) =", mp.nstr(mn,10), " (>=0 needed)")
print("final cell true pairs: min Btrue =", mp.nstr(mn2,10), "; max Btrue =", mp.nstr(mxB,10))
print("code breg_abs ~ 0.25 (box) >= max Btrue:", 0.25 > float(mxB))
# exact symbolic: Btrue = -sig/rho + rho/(2(2+rt2 t)) + (4 sig + 2 y^2 sig^2)/(2 rho (2+rt2 t)) - rho A1'/(y^2 (2+rt2 t))
# verify my exact regular form
mn3 = mp.mpf('9'); mx3 = mp.mpf('-9')
for _ in range(200):
    y = mp.mpf(10)**random.uniform(-9,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    rho = 1+y*y*random.random(); sig=(rho-1)/y**2; x0=1/(rho*y); t=rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t); Ap = 1-mp.atan(zz)/zz
    Btrue = (phi_true + pi/4)/((pi/sqrt2)*y)
    Bmine = -sig/rho + rho/(2*(2+sqrt2*t)) + (4*sig+2*y*y*sig*sig)/(2*rho*(2+sqrt2*t)) - rho*Ap/(y*y*(2+sqrt2*t))
    d = Btrue-Bmine
    mn3 = min(mn3,d); mx3 = max(mx3,d)
print("my exact regular form: min/max residual =", mp.nstr(mn3,8), mp.nstr(mx3,8))
# so D13h breg_h is not exact; but is breg_h - Btrue >= 0 PROVABLY on true saddles?
# breg_h - Btrue = [-sig/(2rho)+sig/rho] + [rho/(2(2+rt2t)) - rho/(2(2+rt2t))]
#   + [(-rt2 y sig + 2 sig)/(2(2+rt2 t)) - (4 sig + 2 y^2 sig^2)/(2 rho(2+rt2 t))]
#   + [-rho A1'/(2+rt2 t) + rho A1'/(y^2 (2+rt2 t))]
print()
print("decompose breg_h - Btrue symbolically:")
import sympy as sp
y_, s_ = sp.symbols('y sigma', positive=True)
r_ = 1+y_**2*s_; t_ = r_*y_
A1_ = sp.symbols('Ap', positive=True)
bh = -s_/(2*r_) + (-2*y_*s_ + 2*sp.sqrt(2)*s_ + sp.sqrt(2)*r_ - 2*sp.sqrt(2)*r_*A1_)/(2*sp.sqrt(2)*(2+sp.sqrt(2)*t_))
Bm = -s_/r_ + r_/(2*(2+sp.sqrt(2)*t_)) + (4*s_+2*y_**2*s_**2)/(2*r_*(2+sp.sqrt(2)*t_)) - r_*A1_/(y_**2*(2+sp.sqrt(2)*t_))
diff = sp.simplify(sp.together(bh - Bm))
print("breg_h - Btrue =", sp.nsimplify(sp.simplify(diff.subs(A1_,0).subs(s_, sp.Rational(1,2)).subs(y_,0))))
num = sp.fraction(sp.together(diff))[0]
num_exp = sp.expand(num)
print("numerator (before A1' substitution):", num_exp)
