"""Certified bounds for |c0|, |c1| (squircle W4 residual constants).
Method: rigorous ball arithmetic (python-flint arb/acb).
  U, V: acb.integral on three holomorphic pieces (subs v=b+e^s, v=a-e^s),
        plus endpoint tails bounded by exact ∫(log(1/s)+C)^k formulas.
  R0, R1: all atoms of conjecture_terms() as balls (z.polylog(s) = Li_s(z)).
Run: python w4_certified.py            # PREC=1100, EPS=2^-1300
"""
import sys, time
sys.dont_write_bytecode = True
from math import comb, factorial
from flint import arb, acb, ctx
from w4_symbol import conjecture_terms

PREC, EPS_EXP = 1100, 1300          # tail < 2^-1267 ~= 1e-381
ctx.prec = PREC
TARGET = arb(10) ** (-300)

def mag(z):                          # ball containing |z'| for z' in z (arb)
    return abs(z.real) + abs(z.imag)
def mul(*xs):                        # avoid arb_pow NaN on zero-centered balls
    out = xs[0]
    for x in xs[1:]:
        out = out * x
    return out

eps   = arb(2) ** (-EPS_EXP)
sqrt2 = acb(2).sqrt()
a, b  = 1 + sqrt2, sqrt2 - 1
log2  = acb(2).log()
pi    = acb.pi()
c1_6, c1_3   = acb(arb(1)/arb(6)), acb(arb(1)/arb(3))
c7_48, c1_16 = acb(arb(7)/arb(48)), acb(arb(1)/arb(16))
c7_48r, c1_16r = arb(7)/arb(48), arb(1)/arb(16)

def fU_mid(v, analytic=False):
    l = ((v - b)*(v + a)/(4*v)).log()
    return l*l*l * (c1_6/v - c1_3/(v + a))
def fV_mid(v, analytic=False):
    l = ((v - b)*(v + a)/(4*v)).log(); m = ((a - v)*(v + b)/(4*v)).log()
    return -c7_48*l*l*l/v + c1_16*l*l*m/v
def fU_left(s, analytic=False):      # v = b + e^s in [b+eps, 1]
    e = s.exp(); v = b + e
    l = s + ((v + a)/(4*v)).log()
    return l*l*l * (c1_6/v - c1_3/(v + a)) * e
def fV_left(s, analytic=False):
    e = s.exp(); v = b + e
    l = s + ((v + a)/(4*v)).log(); m = ((a - v)*(v + b)/(4*v)).log()
    return (-c7_48*l*l*l/v + c1_16*l*l*m/v) * e
def fU_right(s, analytic=False):     # v = a - e^s in [a-1, a-eps]
    e = s.exp(); v = a - e
    l = ((v - b)*(v + a)/(4*v)).log()
    return l*l*l * (c1_6/v - c1_3/(v + a)) * e
def fV_right(s, analytic=False):
    e = s.exp(); v = a - e
    l = ((v - b)*(v + a)/(4*v)).log(); m = s + ((v + b)/(4*v)).log()
    return (-c7_48*l*l*l/v + c1_16*l*l*m/v) * e

s_lo, s_ml, s_mr = -EPS_EXP*log2, (1 - b).log(), acb(0)
t0 = time.time()
IU = (acb.integral(fU_left, s_lo, s_ml) + acb.integral(fU_mid, acb(1), a - 1)
      + acb.integral(fU_right, s_lo, s_mr))
IV = (acb.integral(fV_left, s_lo, s_ml) + acb.integral(fV_mid, acb(1), a - 1)
      + acb.integral(fV_right, s_lo, s_mr))
print('acb.integral time: %.1fs' % (time.time() - t0), flush=True)

# ---- endpoint tails (rigorous sup bounds over covering balls) ----
w = EPS_EXP * arb(2).log()                       # w = log(1/eps)
def P(C, k):                                     # ∫_0^eps (log(1/s)+C)^k ds / eps
    return sum(comb(k, j)*mul(*[w + C]*j)*factorial(k - j) for j in range(k + 1))
vb = acb(b) + acb(arb(0, eps))                   # covers [b, b+eps]
C1, invb = mag(((vb + a)/(4*vb)).log()), mag(1/vb)
H_U_b = mag(c1_6/vb - c1_3/(vb + a))
M1    = mag(((a - vb)*(vb + b)/(4*vb)).log())
T_U_b = H_U_b * eps * P(C1, 3)
T_V_b = c7_48r*invb*eps*P(C1, 3) + c1_16r*invb*M1*eps*P(C1, 2)
va = acb(a) + acb(arb(0, eps))                   # covers [a-eps, a]
C2, inva = mag(((va + b)/(4*va)).log()), mag(1/va)
L1    = mag(((va - b)*(va + a)/(4*va)).log())
H_U_a = mag(c1_6/va - c1_3/(va + a))
S1    = (w + C2) + 1
T_U_a = H_U_a * eps * mul(L1, L1, L1)
T_V_a = c7_48r*inva*eps*mul(L1, L1, L1) + c1_16r*inva*eps*mul(L1, L1)*S1
rU, rV = T_U_b + T_U_a, T_V_b + T_V_a
U = IU + acb(arb(0, rU))
V = IV + acb(arb(0, rV))
print('U =', U, '\n  (integral rad %.3e, tail rad %.3e)' % (mag(IU-U+U-IU), rU))
print('V =', V)

# ---- R0, R1 from the literal conjecture terms ----
def qs2(p, q):                                   # p + q*sqrt(2), Fractions
    v = acb(arb(p.numerator)/arb(p.denominator))
    return v + acb(arb(q.numerator)/arb(q.denominator))*sqrt2 if q else v
def atom_value(at):
    if at[0] == 'log':
        z = qs2(at[1].p, at[1].q)
        assert z.real > 0
        return acb(abs(z.real)).log()
    if at[0] == 'li':
        return qs2(at[2].p, at[2].q).polylog(acb(at[1]))
    if at == ('pi2',):
        return pi*pi
    assert at == ('z3',)
    return acb(3).zeta()
R0, R1 = acb(0), acb(0)
for c, atoms in conjecture_terms():
    val = acb(1)
    for at in atoms:
        val = val * atom_value(at)
    R0 += acb(arb(c.p.numerator)/arb(c.p.denominator)) * val
    R1 += acb(arb(c.q.numerator)/arb(c.q.denominator)) * val

# ---- residuals and the certified test ----
pi4 = pi*pi*pi*pi
c0, c1 = (U - R0)/pi4, (V - R1)/pi4
for name, ci in [('c0', c0), ('c1', c1)]:
    print(name, '=', ci)
    print('  |%s| bound ball:' % name, mag(ci))
    ok = bool(mag(ci) < TARGET)     # certified ball comparison
    print('  CERTIFIED |%s| < 1e-300 :' % name, ok)
    assert ok

# ---- NUMERICAL cross-check (mpmath, not part of the certificate) ----
import mpmath as mp
mp.mp.dps = 100
def num(x): return mp.mpf(str(x.p)) + mp.sqrt(2)*mp.mpf(str(x.q))
am, bm = num(__import__('w4_symbol').A), num(__import__('w4_symbol').B)
fU = lambda v: mp.log((v-bm)*(v+am)/(4*v))**3*(1/(6*v)-1/(3*(v+am)))
fV = lambda v: (-7*mp.log((v-bm)*(v+am)/(4*v))**3/(48*v)
                + mp.log((v-bm)*(v+am)/(4*v))**2*mp.log((am-v)*(v+bm)/(4*v))/(16*v))
print('NUMERICAL mp  U =', mp.nstr(mp.quad(fU, [bm, 1, am-1, am]), 50))
print('NUMERICAL mp  V =', mp.nstr(mp.quad(fV, [bm, 1, am-1, am]), 50))
