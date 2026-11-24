"""Referee step 1: M_2, M_3 from the definition, independently of the project's scripts.

Definition (main.tex eq. F, Lemma holo): P/4 = 2 - eps F(eps), F = int_0^1 Phi, M_k = -[eps^k]F / 2.
On 0<t<=1/2 we use the cancellation-free form from the proof of Lemma holo:
  Phi = B e^{eps y} - B^2 e^{2 eps y} / (A e^{eps x} (1 + sqrt(1+z))),  z = (B e^{eps y}/(A e^{eps x}))^2,
  A = 1/t, B = 1/(1-t), x = log t, y = log(1-t), and F = 2 int_0^{1/2} Phi.
Usage: python s1_numeric.py DPS
"""
import sys
import sympy as sp
import mpmath as mp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 90
mp.mp.dps = DPS

e, r, x, y, B = sp.symbols('e r x y B', positive=True)
# Phi / B written with r = B/A = t/(1-t)
z = r**2 * sp.exp(2*e*(y - x))
phiB = sp.exp(e*y) - r*sp.exp(e*(2*y - x)) / (1 + sp.sqrt(1 + z))
ser = sp.series(phiB, e, 0, 4).removeO()
coef = {k: sp.simplify(ser.coeff(e, k)) for k in range(4)}
fk = {k: sp.lambdify((r, x, y), coef[k], 'mpmath') for k in range(4)}


def Mk(k):
    def integrand(t):
        t = mp.mpf(t)
        return fk[k](t/(1-t), mp.log(t), mp.log(1-t)) / (1-t)
    F_k = 2*mp.quad(integrand, [0, mp.mpf(1)/4, mp.mpf(1)/2])
    return -F_k/2


s2 = mp.sqrt(2)
a = 1 + s2
b = s2 - 1
L = mp.log(a)
lam = mp.log(2)
pi2 = mp.pi**2
z3 = mp.zeta(3)
Li3 = lambda u: mp.polylog(3, u)
Li2 = lambda u: mp.polylog(2, u)

T = (pi2 + 6*lam**2)/24

# (T1) transcribed from main.tex
M2_T1 = (T + (189*s2/64 - mp.mpf(19)/16)*z3 + (11 + 15*s2)*Li3(b) + (11 - 15*s2)*Li3(-b)
         - 27*s2*Li3(1/s2) + (12*L - 7*s2/2*lam)*Li2(b)
         + (15*s2/8 - mp.mpf(7)/12)*pi2*L - (mp.mpf(1)/24 + 65*s2/48)*pi2*lam
         + (mp.mpf(1)/6 + 9*s2/16)*lam**3 + 11*s2/8*lam**2*L - (1 + 7*s2/4)*lam*L**2
         + (mp.mpf(25)/6 - 7*s2/3)*L**3)

# (T2) transcribed from main.tex
W3p_T2 = ((7*s2/32 - mp.mpf(3)/16)*z3 + (1 + s2)*Li3(b) + (1 - s2)*Li3(-b) - 2*s2*Li3(1/s2)
          + (2*L - s2*lam)*Li2(b) + (s2/12 - mp.mpf(1)/8)*pi2*L - (1 + s2)/24*pi2*lam
          + (1 + s2)/24*lam**3 + s2/4*lam**2*L - (1 + s2)/2*lam*L**2 + (5 - s2)/6*L**3)


def l1(v):
    return mp.log((v - b)*(v + a)/(4*v))


def l2(v):
    return mp.log((a - v)*(v + b)/(4*v))


def vint(f):
    return mp.quad(f, [b, 1, a])


# W_4 as printed in thm:M3
W4 = vint(lambda v: l1(v)**3*((8 - 7*s2)/(48*v) - 1/(3*(v + a))) + s2/16*l1(v)**2*l2(v)/v)
# T as an integral (app:M2)
T_int = vint(lambda v: l1(v)*(1/(2*v) + 1/(2*(v - a)) - 1/(v + a) - 1/(2*(v + b))))
# W_3 and W_3' as integrals (squircle note eq. W3, W4 block) -- used only as cross-checks
W3_int = vint(lambda v: l1(v)**2*((4 - 3*s2)/(8*v) - 1/(v + a)) + s2/8*l1(v)*l2(v)/v)
W3p_int = vint(lambda v: l1(v)**2*((1 - s2)/(4*v) + 1/(4*(v - a)) - 1/(4*(v + a)) - 1/(2*(v + b)))
               + l1(v)*l2(v)*((1 + s2)/(4*v) - 1/(2*(v + a))))

M1 = Mk(1)
M2 = Mk(2)
M3 = Mk(3)
M1_thm = ((19*s2/48 - mp.mpf(1)/24)*pi2 + lam**2/2 + 2*s2*lam*L - (1 + 5*s2/2)*L**2 - 5*s2*Li2(b))


def dig(u, v):
    d = abs(u - v)
    return mp.inf if d == 0 else -mp.log10(d)


nd = 40
print('dps', DPS)
print('sanity M1 def      ', mp.nstr(M1, nd), ' vs thm:M1 digits', mp.nstr(dig(M1, M1_thm), 5))
print('M2 (definition)    ', mp.nstr(M2, nd))
print('M2 (T1 closed form)', mp.nstr(M2_T1, nd), ' digits', mp.nstr(dig(M2, M2_T1), 5))
print('M2 printed -0.02216191108726939942 digits',
      mp.nstr(dig(M2, mp.mpf('-0.02216191108726939942')), 5))
print('T integral vs (pi^2+6lam^2)/24 digits', mp.nstr(dig(T_int, T), 5))
print('M2 - T - W3_int digits', mp.nstr(dig(M2, T + W3_int), 5))
print('M3 (definition)    ', mp.nstr(M3, nd))
print('W4 (printed integral)', mp.nstr(W4, nd))
print('W3p (T2)           ', mp.nstr(W3p_T2, nd), ' vs W3p integral digits', mp.nstr(dig(W3p_T2, W3p_int), 5))
print('M3 vs W4+W3p(T2)+T digits', mp.nstr(dig(M3, W4 + W3p_T2 + T), 5))
print('M3 printed 0.02277052149178657746 digits', mp.nstr(dig(M3, mp.mpf('0.02277052149178657746')), 5))
print('W4 printed 0.79578505983224081132 digits', mp.nstr(dig(W4, mp.mpf('0.79578505983224081132')), 5))
# side identities quoted in the text around thm:M2
chi3 = (Li3(b) - Li3(-b))/2
Lchi8 = (mp.zeta(3, mp.mpf(1)/8) - mp.zeta(3, mp.mpf(3)/8) - mp.zeta(3, mp.mpf(5)/8) + mp.zeta(3, mp.mpf(7)/8))/512
print('L(chi8,3)', mp.nstr(Lchi8, 15))
print('eq chi3 digits', mp.nstr(dig(48*chi3, 32*s2*Lchi8 - 3*pi2*L + 4*L**3), 5))
print('Li3(1/2) formula digits', mp.nstr(dig(Li3(mp.mpf(1)/2), mp.mpf(7)/8*z3 - pi2/12*lam + lam**3/6), 5))
print('T = (zeta2-Li2(1/2))/2 digits', mp.nstr(dig(T, (pi2/6 - Li2(mp.mpf(1)/2))/2), 5))
