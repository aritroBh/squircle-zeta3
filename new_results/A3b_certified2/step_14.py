
script = r'''"""Certified bounds for |c0|, |c1| (squircle W4 residual constants).
Method: rigorous ball arithmetic (python-flint arb/acb).
  U, V: acb.integral on three holomorphic pieces (subs v=b+e^s, v=a-e^s),
        plus endpoint tails bounded by exact formulas for int(log(1/s)+C)^k ds.
  R0, R1: all atoms of conjecture_terms() as balls (z.polylog(s) = Li_s(z)).
RIGOR (analytic flag): acb.integral calls f(x, analytic); with analytic=True
  the integrand must be non-finite whenever it is not holomorphic on the ball,
  else the Gauss-Legendre error bound is not rigorous. Every branch-cut
  function in the six integrands is log; each is called as .log(analytic=analytic).
  (exp is entire; division is meromorphic - poles give non-finite automatically.)
Run: python w4_certified.py [PREC=1100] [EPS_EXP=PREC+100]
"""
import sys, time
sys.dont_write_bytecode = True
from math import comb, factorial
from flint import arb, acb, ctx
from w4_symbol import conjecture_terms

PREC = int(sys.argv[1]) if len(sys.argv) > 1 else 1100
EPS_EXP = int(sys.argv[2]) if len(sys.argv) > 2 else PREC + 100
ctx.prec = PREC
TARGET = arb(10) ** (-300)
print('PREC =', PREC, ' EPS_EXP =', EPS_EXP, flush=True)

def mag(z):                          # arb ball containing |z'| for z' in z
    return abs(z.real) + abs(z.imag)

eps   = arb(2) ** (-EPS_EXP)
sqrt2 = acb(2).sqrt()
a, b  = 1 + sqrt2, sqrt2 - 1
log2  = acb(2).log()
pi    = acb.pi()
c1_6, c1_3   = acb(arb(1)/arb(6)), acb(arb(1)/arb(3))
c7_48, c1_16 = acb(arb(7)/arb(48)), acb(arb(1)/arb(16))
c7_48r, c1_16r = arb(7)/arb(48), arb(1)/arb(16)

# ---- integrands: analytic flag forwarded to EVERY branch-cut call (log) ----
def fU_mid(v, analytic=False):
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    return l*l*l * (c1_6/v - c1_3/(v + a))
def fV_mid(v, analytic=False):
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    m = ((a - v)*(v + b)/(4*v)).log(analytic=analytic)
    return -c7_48*l*l*l/v + c1_16*l*l*m/v
def fU_left(s, analytic=False):      # v = b + e^s in [b+eps, 1]
    e = s.exp(); v = b + e
    l = s + ((v + a)/(4*v)).log(analytic=analytic)
    return l*l*l * (c1_6/v - c1_3/(v + a)) * e
def fV_left(s, analytic=False):
    e = s.exp(); v = b + e
    l = s + ((v + a)/(4*v)).log(analytic=analytic)
    m = ((a - v)*(v + b)/(4*v)).log(analytic=analytic)
    return (-c7_48*l*l*l/v + c1_16*l*l*m/v) * e
def fU_right(s, analytic=False):     # v = a - e^s in [a-1, a-eps]
    e = s.exp(); v = a - e
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    return l*l*l * (c1_6/v - c1_3/(v + a)) * e
def fV_right(s, analytic=False):
    e = s.exp(); v = a - e
    l = ((v - b)*(v + a)/(4*v)).log(analytic=analytic)
    m = s + ((v + b)/(4*v)).log(analytic=analytic)
    return (-c7_48*l*l*l/v + c1_16*l*l*m/v) * e
# branch-cut audit: integrands contain only log (flagged), exp (entire),
# +,-,*,/ (meromorphic). No sqrt, no pow anywhere (powers are products).

s_lo, s_ml, s_mr = -EPS_EXP*log2, (1 - b).log(), acb(0)
t0 = time.time()
IU_left  = acb.integral(fU_left,  s_lo, s_ml)
IU_mid   = acb.integral(fU_mid,   acb(1), a - 1)
IU_right = acb.integral(fU_right, s_lo, s_mr)
IV_left  = acb.integral(fV_left,  s_lo, s_ml)
IV_mid   = acb.integral(fV_mid,   acb(1), a - 1)
IV_right = acb.integral(fV_right, s_lo, s_mr)
IU = IU_left + IU_mid + IU_right
IV = IV_left + IV_mid + IV_right
print('acb.integral time: %.1fs' % (time.time() - t0), flush=True)
print('IU_left  =', IU_left,  '\nIU_mid   =', IU_mid,  '\nIU_right =', IU_right)
print('IV_left  =', IV_left,  '\nIV_mid   =', IV_mid,  '\nIV_right =', IV_right)
for name, I in [('IU', IU), ('IV', IV)]:
    print('%s = %s\n  integral rad: re %.3e im %.3e' % (name, I, I.real.rad(), I.imag.rad()))
assert IU.is_finite() and IV.is_finite(), 'integral not finite'

# ---- endpoint tails (rigorous sup bounds over covering balls) ----
w = EPS_EXP * arb(2).log()                       # w = log(1/eps)
def P(C, k):                                     # (1/eps) * int_0^eps (log(1/s)+C)^k ds
    out = arb(0)                                 #   = sum_j C(k,j) (w+C)^j (k-j)!
    for j in range(k + 1):
        t = arb(comb(k, j)) * arb(factorial(k - j))
        for _ in range(j):                       # x*x*x, not x**j (arb_pow NaN on 0-balls)
            t = t * (w + C)
        out = out + t
    return out
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
T_U_a = H_U_a * eps * L1*L1*L1
T_V_a = c7_48r*inva*eps*L1*L1*L1 + c1_16r*inva*eps*L1*L1*S1
print('tail factors: C1=%s invb=%s H_U_b=%s M1=%s' % (C1, invb, H_U_b, M1))
print('tail factors: C2=%s inva=%s L1=%s H_U_a=%s S1=%s' % (C2, inva, L1, H_U_a, S1))
print('T_U_b =', T_U_b, '\nT_V_b =', T_V_b, '\nT_U_a =', T_U_a, '\nT_V_a =', T_V_a)
rU, rV = T_U_b + T_U_a, T_V_b + T_V_a
print('tail rad rU =', rU, '\ntail rad rV =', rV)
U = IU + acb(arb(0, rU.abs_upper()))
V = IV + acb(arb(0, rV.abs_upper()))
print('U =', U, '\n  tail rad %.3e' % rU.abs_upper())
print('V =', V, '\n  tail rad %.3e' % rV.abs_upper())

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
print('R0 =', R0, '\n  rad: re %.3e im %.3e' % (R0.real.rad(), R0.imag.rad()))
print('R1 =', R1, '\n  rad: re %.3e im %.3e' % (R1.real.rad(), R1.imag.rad()))

# ---- residuals and the certified test ----
pi4 = pi*pi*pi*pi
c0, c1 = (U - R0)/pi4, (V - R1)/pi4
def certified_decimal_exp(ci):
    """Maximal integer k with mag(ci) < 10^-k certifiably (arb < semantics)."""
    m = mag(ci)
    t = -(m.abs_upper().log() / arb(10).log())   # ball approximating -log10|ci|
    k = int(t.mid().floor())
    while not bool(m < arb(10) ** (-k)):
        k -= 1
    while bool(m < arb(10) ** (-(k + 1))):
        k += 1
    return k
for name, ci in [('c0', c0), ('c1', c1)]:
    print(name, '=', ci)
    print('  |%s| bound ball:' % name, mag(ci))
    k = certified_decimal_exp(ci)
    okk = bool(mag(ci) < arb(10) ** (-k))
    ok3 = bool(mag(ci) < TARGET)                # certified ball comparison
    print('  CERTIFIED |%s| < 10^-%d :' % (name, k), okk)
    print('  CERTIFIED |%s| < 1e-300 :' % name, ok3)
    assert okk
    if PREC >= 1100:
        assert ok3, 'production run must beat 1e-300'

# ---- NUMERICAL cross-check (mpmath, not part of the certificate) ----
import mpmath as mp
mp.mp.dps = 100
def fmp(fr): return mp.mpf(fr.numerator) / mp.mpf(fr.denominator)
def num(x): return fmp(x.p) + mp.sqrt(2)*fmp(x.q)
am, bm = num(__import__('w4_symbol').A), num(__import__('w4_symbol').B)
fU = lambda v: mp.log((v-bm)*(v+am)/(4*v))**3*(1/(6*v)-1/(3*(v+am)))
fV = lambda v: (-7*mp.log((v-bm)*(v+am)/(4*v))**3/(48*v)
                + mp.log((v-bm)*(v+am)/(4*v))**2*mp.log((am-v)*(v+bm)/(4*v))/(16*v))
print('NUMERICAL mp  U =', mp.nstr(mp.quad(fU, [bm, 1, am-1, am]), 50))
print('NUMERICAL mp  V =', mp.nstr(mp.quad(fV, [bm, 1, am-1, am]), 50))
'''
with open('w4_certified.py', 'w') as f:
    f.write(script)
print("written", len(script), "bytes")
import subprocess, sys
print(subprocess.run([sys.executable, '-c', 'import ast; ast.parse(open("w4_certified.py").read()); print("syntax OK")'],
                     capture_output=True, text=True).stdout)
