# step_10.py : Verify Mellin identity and recompute A_D, B_D, C_D by quadrature at high precision
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta, gamma, sinh, cosh, exp, tanh, inf, quad
mp.dps = 100
L2 = log(2); t1 = log(1+sqrt(2))

def h(w):
    w = mpf(w)
    return -1 - w/2 + (2+w*w)/(2*sqrt(1+w*w))
def q(w):
    return h(w)/w

# --- quadrature definitions (same as prev) ---
A_D = quad(lambda w: q(w)+mpf('0.5'), [0,1]) + quad(lambda w: q(w)+1/w, [1, inf])
B_D = quad(lambda w: (q(w)+mpf('0.5'))*log(w), [0,1]) + quad(lambda w: (q(w)+1/w)*log(w), [1, inf])
C_D = quad(lambda w: (q(w)+mpf('0.5'))*log(w)**2, [0,1]) + quad(lambda w: (q(w)+1/w)*log(w)**2, [1, inf])
print("A_D quad =", nstr(A_D, 60))
print("B_D quad =", nstr(B_D, 60))
print("C_D quad =", nstr(C_D, 60))

# --- my closed-form candidates from psi-expansion ---
z2 = zeta(2); z3 = zeta(3)
A_cf = L2
B_cf = L2**2/2 + z2/2 - L2/2
C_cf = L2**3/3 + L2*z2 + z3/2 - L2**2/2 - z2/2 + L2
print("A_cf =", nstr(A_cf, 60), " diff:", nstr(A_D-A_cf, 6))
print("B_cf =", nstr(B_cf, 60), " diff:", nstr(B_D-B_cf, 6))
print("C_cf =", nstr(C_cf, 60), " diff:", nstr(C_D-C_cf, 25))

# --- numerical check of Mellin identity M_h(z) = G1(z)(2+z)/(1+z) in the strip ---
def Mh_num(z):
    # ∫_0^∞ w^{z-1} h(w) dw, valid -1<Re z<0
    z = mpf(z)
    return quad(lambda w: w**(z-1)*h(w), [0, 1, inf])
def G1(z):
    z = mpf(z)
    return gamma(z/2)*gamma(mpf('0.5')-z/2)/(4*sqrt(pi))
for z in ['-0.5','-0.3','-0.7','-0.25']:
    lhs = Mh_num(z); rhs = G1(z)*(2+mpf(z))/(1+mpf(z))
    print("z=%s  Mh_num=%s  G1(2+z)/(1+z)=%s  diff=%s" % (z, nstr(lhs,25), nstr(rhs,25), nstr(lhs-rhs,6)))
