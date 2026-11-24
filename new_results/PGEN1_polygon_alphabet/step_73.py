# step_45: decagon canonical pair: N3 = a*zeta(3) + b*L(chi_5,3), via single-valued L3
import mpmath as mp
mp.mp.dps = 80
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    """Zagier single-valued trilogarithm, real x < 1:
    L3(x) = Li3(x) - log|x| Li2(x) + (1/3) log^2|x| Li1(x),  Li1(x) = -log(1-x)."""
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
def L3r(x):
    """L3 for any real x != 0,1: use inversion L3(x) = L3(1/x)."""
    x = mp.mpf(x)
    if x < 1: return L3(x)
    return L3(1/x)
# controls for L3 functional equations
z3 = mp.zeta(3)
print("L3(1) =", L3(1) if False else "n/a; zeta3 =", z3)
print("L3(1/2) =", L3(mp.mpf('0.5')), " 7/8 zeta3 =", mp.mpf(7)/8*z3, " diff:", L3(mp.mpf('0.5'))-mp.mpf(7)/8*z3)
# Landen: L3(x)+L3(1-x)+L3(1-1/x) = zeta(3)
for xx in ['0.3','0.7','-0.4','2.5']:
    x = mp.mpf(xx)
    s = L3r(x)+L3r(1-x)+L3r(1-1/x)
    print(f"  Landen L3 at {xx}: {s}  diff {s-z3}")
# duplication: L3(x^2) = 4L3(x)+4L3(-x)
for xx in ['0.6','-0.8']:
    x = mp.mpf(xx)
    print(f"  dupl L3 at {xx}:", L3r(x**2)-4*L3r(x)-4*L3r(-x))
# L(chi5,3)
phi = (1+mp.sqrt(5))/2
L5 = mp.nsum(lambda n: mp.kronecker(5,n)/mp.mpf(n)**3, [1, mp.inf])
print("L(chi5,3) =", L5)
X1, Xm1 = L3(1/phi), L3(-1/phi)
alpha = (X1+Xm1)/(2*z3)
beta  = (X1-Xm1)/(2*L5)
print("alpha =", alpha, "  beta =", beta)
import sympy as sp
print("alpha nsimplify:", sp.nsimplify(str(alpha)), " beta nsimplify:", sp.nsimplify(str(beta)))
