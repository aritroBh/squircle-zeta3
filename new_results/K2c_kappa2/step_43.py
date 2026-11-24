work='<project>/research/08_theory/K2c_kappa2/'
code = r'''
import sympy as sp
sig = sp.Symbol('sig'); q = sp.Symbol('q'); Z = q**2; g1 = sp.Symbol('g1')
s2 = sp.sqrt(2); ORD = 6
def T(e, n=ORD): return sp.series(e, sig, 0, n).removeO()
def saddle_p(p):
    sp_ = sp.sqrt(p)*sig
    yp = T(sp.expand(T(sp.sqrt(1 + q**2*sp_**2/4)) - q*sp_/2))
    return sp_, yp, T(sp.expand(q*yp/sp_))
def phi_p(p):
    sp_, yp, xp = saddle_p(p)
    return T(sp.expand(p*q*yp/sp_ - 2*sp.log(1+(yp-1))/sig**2))
def phi2_p(p):
    sp_, yp, xp = saddle_p(p)
    return T(sp.expand(-sig**-2*((xp+Z)**(-2) - xp**(-2))))
lead = 2*(s2-1)*q/sig
D = T(sp.expand(phi_p(2) - phi_p(1)))
S_exp = T(sp.exp(T(sp.expand(-(D - lead)))))
sp1, y1, x1 = saddle_p(1); sp2, y2, x2 = saddle_p(2)
w1 = T(sp.expand(x1/y1**2)); w2 = T(sp.expand(x2/y2**2))
Pref = T(sp.expand((phi2_p(1)/phi2_p(2))**(sp.Rational(1,2))*(x1/x2)*(w2/w1)**(sp.Rational(3,2))))
# CORRECTED bracket: effective parameter sig/sqrt(p)
def B0(s): return 1 + s/(2*q) - sp.Rational(5,8)*s**2
def B1(s): return 1 + 3*s/(2*q) + (-sp.Rational(7,4) + sp.Rational(3,4)/q**2)*s**2
Br = T(sp.expand(B0(sig/s2)/B0(sig)))
Nf = sig**-2
m1_1 = T(sp.expand(sp.Rational(3,2)*g1*w1/Nf * B1(sig)/B0(sig)))
m1_2 = T(sp.expand(sp.Rational(3,2)*g1*w2/Nf * B1(sig/s2)/B0(sig/s2)))
F4 = T(sp.expand((1+m1_2)/(1+m1_1)))
tot = T(sp.expand(sp.Rational(-2)*S_exp*Pref*Br*F4))
totn = T(sp.expand(tot/sp.exp(Z/2)))
d1 = sp.expand(totn).coeff(sig, 1)
print("R*e^lead/e^{Z/2} =", sp.expand(totn))
print("d1 =", sp.simplify(d1))
# numeric at the two angles vs quadrature fit
import mpmath
mpmath.mp.dps = 60
d1f = sp.lambdify([q, g1], d1, 'mpmath')
for av, gv in [(mpmath.pi/2, -mpmath.mpf('0.5')), (5*mpmath.pi/6, None), (mpmath.pi/3, None), (2*mpmath.pi/3, None), (3*mpmath.pi/4, None)]:
    a = mpmath.mpf(av)
    if gv is None:
        gv = -mpmath.mpf(3)/4 - 1j*mpmath.e**(1j*a)/(4*mpmath.sin(a))
    qv = mpmath.sqrt(1j*a)
    print('alpha=%s: d1 = %s' % (mpmath.nstr(a,5), mpmath.nstr(d1f(qv, gv), 18)))
'''
open(work+'step_32.py','w').write(code)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', work+'step_32.py'],capture_output=True,text=True,timeout=900)
print(r.stdout[-3000:]); print(r.stderr[-500:])
