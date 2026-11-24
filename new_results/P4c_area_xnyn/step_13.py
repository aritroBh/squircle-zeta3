# step_11.py : high-precision A_D,B_D,C_D via t = asinh(w) (well-conditioned) + closed forms
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta, sinh, cosh, exp, tanh, inf, quad
mp.dps = 120
L2 = log(2); t1 = log(1+sqrt(2)); z2 = zeta(2); z3 = zeta(3)

def F0(t): return tanh(t/2)*(cosh(t)-1)/2        # = (q+1/2) dw/dt , t in [0,t1]
def F1(t): return (3+exp(-2*t))/(4*sinh(t))       # = (q+1/w) dw/dt , t in [t1,inf)

def t_int(weight):
    I0 = quad(lambda t: F0(t)*weight(log(sinh(t))), [0, t1/2, t1])
    I1 = quad(lambda t: F1(t)*weight(log(sinh(t))), [t1, 2*t1, t1+10, t1+40, t1+120])
    return I0 + I1
# tail bound beyond t1+120: |F1| <= (3+e^{-2t})/(4 sinh t) <= e^{-t}, weight <= (t)^2 -> negligible at 120 dps
A_D = t_int(lambda u: 1)
B_D = t_int(lambda u: u)
C_D = t_int(lambda u: u*u)
A_cf = L2
B_cf = L2**2/2 + z2/2 - L2/2
C_cf = L2**3/3 + L2*z2 + z3/2 - L2**2/2 - z2/2 + L2
print("A_D  =", nstr(A_D, 70)); print("  closed =", nstr(A_cf, 70)); print("  diff =", nstr(A_D-A_cf, 8))
print("B_D  =", nstr(B_D, 70)); print("  closed =", nstr(B_cf, 70)); print("  diff =", nstr(B_D-B_cf, 8))
print("C_D  =", nstr(C_D, 70)); print("  closed =", nstr(C_cf, 70)); print("  diff =", nstr(C_D-C_cf, 8))

# Assemble R3(n) = (1/3) log^3 n + b3 log^2 n + c3 log n + d3  (prev's recursion, L = log2/2)
L = L2/2
b3 = L + mpf('0.5') - A_cf
c3 = L**2 - 2*L*(A_cf - mpf('0.5')) + 1 + 2*B_cf
d3 = L**3/3 - L**2*(A_cf - mpf('0.5')) + L*(1+2*B_cf) + 1 - C_cf + 5*sqrt(2)/4 + 1/(2*sqrt(2))
print("b3 =", nstr(b3, 60)); print("c3 =", nstr(c3, 60)); print("d3 =", nstr(d3, 60))
# simplified closed forms
b3s = mpf('0.5') - L2/2
c3s = L2**2/4 - L2/2 + 1 + z2
d3s = -L2**3/24 + L2**2/8 - L2/2 - L2*z2/2 + z2/2 - z3/2 + 1 + 3*sqrt(2)/2
print("b3 simplified diff:", nstr(b3-b3s, 8))
print("c3 simplified diff:", nstr(c3-c3s, 8))
print("d3 simplified diff:", nstr(d3-d3s, 8))
