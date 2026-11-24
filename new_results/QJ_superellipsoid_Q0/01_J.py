# 01_J.py -- exact derivation of the closed form of
#   J(q) = int_0^oo [1/(1+r+R) - 1/(1+q)] dr/(r(1+r)),   R = sqrt(q^2+r^2),  0 < q <= 1,
# by the Euler substitution u = r + R (rational in u), partial fractions, and an antiderivative
# that is verified by exact differentiation.  Numerical check to 70 digits at several q.
import sympy as sp
import mpmath as mp

q, u, r, be = sp.symbols('q u r beta', positive=True)
ok = lambda e: sp.simplify(e) == 0

# --- Euler substitution: r = (u^2-q^2)/(2u), R = (u^2+q^2)/(2u)
rr = (u**2 - q**2)/(2*u)
RR = (u**2 + q**2)/(2*u)
assert ok(RR**2 - (q**2 + rr**2)), "R^2 = q^2 + r^2"
assert ok(rr + RR - u), "r + R = u"
drdu = sp.diff(rr, u)
assert ok(drdu - (u**2 + q**2)/(2*u**2))
print("[J1] r=(u^2-q^2)/(2u), R=(u^2+q^2)/(2u): R^2=q^2+r^2, r+R=u, dr/du=(u^2+q^2)/(2u^2)  (exact)")

integrand_r = (1/(1 + rr + RR) - 1/(1 + q))/(rr*(1 + rr))
integrand_u = sp.together(integrand_r*drdu)
claimed = -2*(u**2 + q**2)/((1 + q)*(1 + u)*(u + q)*(u**2 + 2*u - q**2))
assert ok(integrand_u - claimed)
print("[J2] integrand in u = -2(u^2+q^2)/((1+q)(1+u)(u+q)(u^2+2u-q^2))  (exact)")

# --- factor u^2+2u-q^2 = (u+1-beta)(u+1+beta), beta = sqrt(1+q^2)
B = sp.sqrt(1 + q**2)
assert ok((u + 1 - B)*(u + 1 + B) - (u**2 + 2*u - q**2))
# partial fractions with beta as an independent symbol (valid identity in u for distinct poles),
# then specialise beta = sqrt(1+q^2)
gen = -2*(u**2 + q**2)/((1 + q)*(1 + u)*(u + q)*(u + 1 - be)*(u + 1 + be))
pf = sp.apart(gen, u)
terms = sp.Add.make_args(pf)
poles = {}
for t in terms:
    n, d = sp.fraction(sp.factor(t))
    # d = const * (u - u_k)
    dpoly = sp.Poly(d, u)
    assert dpoly.degree() == 1
    lc = dpoly.LC(); uk = -dpoly.TC()/lc
    poles[sp.simplify(uk)] = sp.simplify(n/lc)
print("[J3] partial fractions: poles", list(poles.keys()))
for k, v in poles.items():
    print("      residue at u =", k, ":", sp.factor(v))
recon = sum(c/(u - k) for k, c in poles.items())
assert ok(recon - gen)
csum = sp.simplify(sum(poles.values()).subs(be, B))
print("[J4] sum of residues (beta=sqrt(1+q^2)) =", csum)
assert csum == 0

# --- antiderivative F(u) = sum c_k log(u - u_k); all u - u_k > 0 for u >= q (checked below)
F = sum(c*sp.log(u - k) for k, c in poles.items())
assert ok(sp.diff(F, u) - gen)
print("[J5] F(u) = sum_k c_k log(u-u_k) ; F'(u) - integrand = 0  (exact, beta symbolic)")
# positivity of u - u_k on [q, oo): u_k in {-1, -q, -1+beta, -1-beta}; q+1-beta > 0 since (q+1)^2-beta^2 = 2q
assert ok((q + 1)**2 - B**2 - 2*q)
print("      u-u_k >= q-u_k > 0 for every pole: q+1>0, 2q>0, q+1-beta=2q/(q+1+beta)>0, q+1+beta>0")
# F(U) -> 0 as U -> oo because sum c_k = 0:  F(U) = sum c_k log(1 - u_k/U) -> 0
# J = F(oo) - F(q) = -F(q).  Work with the logarithms as independent symbols ("letters") so that
# no log-merging by the CAS is involved:  l_k = log(q - u_k).
lg2, lgq, l1, lqb, l1b, lm, lp = sp.symbols('log2 logq log1pq logqpb log1pb logqp1mb logqp1pb')
letter = {-1: l1, -q: lgq + 0, be - 1: lm, -be - 1: lp}
# q - u_k:  q+1, 2q, q+1-beta, q+1+beta
Jder = 0
for k, c in poles.items():
    if k == -q:
        Jder += -c.subs(be, B)*(lg2 + lgq)          # log(2q)
    else:
        Jder += -c.subs(be, B)*letter[k]
print("[J6] J(q) = -F(q) = sum over letters:")
Jder = sp.expand(Jder)
for sym in (l1, lg2, lgq, lm, lp):
    print("      coeff of", sym, ":", sp.factor(sp.radsimp(sp.simplify(Jder.coeff(sym)))))
J1p = sp.expand(Jder.subs(lm, lg2 + lgq - lp))
print("     after log(q+1-beta) = log(2q) - log(q+1+beta):")
for sym in (l1, lg2, lgq, lp):
    print("      coeff of", sym, ":", sp.factor(sp.radsimp(sp.simplify(J1p.coeff(sym)))))
# letter relations (exact algebraic identities between positive quantities, then take log):
assert ok((q + 1 - B)*(q + 1 + B) - 2*q)            # => lm = log2 + logq - lp
assert ok((q + 1 + B)**2 - 2*(q + B)*(1 + B))        # => lp = (log2 + logqpb + log1pb)/2
J1 = sp.expand(Jder.subs(lm, lg2 + lgq - lp))
J1 = sp.expand(J1.subs(lp, (lg2 + lqb + l1b)/2))
# target: [2 log(1+q) - (1+q)(log2+logq)]/(1-q^2) - (log(q+beta) + log(1+beta) - logq)/beta
#   using asinh q = log(q+beta), asinh(1/q) = log((1+beta)/q)
assert ok(sp.asinh(q).rewrite(sp.log) - sp.log(q + B))
assert ok(sp.asinh(1/q).rewrite(sp.log) - sp.log((1 + B)/q))
Jt = (2*l1 - (1 + q)*(lg2 + lgq))/(1 - q**2) - (lqb + l1b - lgq)/B
dd = sp.expand(J1 - Jt)
res = {sym: sp.simplify(dd.coeff(sym)) for sym in (l1, lg2, lgq, lqb, l1b)}
print("[J7] after (q+1-b)(q+1+b)=2q and (q+1+b)^2=2(q+b)(1+b): coefficient differences", res)
assert all(v == 0 for v in res.values()) and sp.simplify(dd.subs({l1: 0, lg2: 0, lgq: 0, lqb: 0, l1b: 0})) == 0
print("     ==> J(q) = [2log(1+q)-(1+q)log(2q)]/(1-q^2) - [asinh q + asinh(1/q)]/sqrt(1+q^2)  PROVED for 0<q<1")

# --- numerical check (70 digits)
mp.mp.dps = 80
def J_num(qq):
    qq = mp.mpf(qq)
    f = lambda t: (1/(1 + t + mp.sqrt(qq*qq + t*t)) - 1/(1 + qq))/(t*(1 + t))
    return mp.quad(f, [0, qq, 1, 10, mp.inf])
def J_cf(qq):
    qq = mp.mpf(qq)
    return (2*mp.log(1 + qq) - (1 + qq)*mp.log(2*qq))/(1 - qq*qq) - (mp.asinh(qq) + mp.asinh(1/qq))/mp.sqrt(1 + qq*qq)
for qq in ['0.7071067811865475244', '0.75', '0.85', '0.95', '0.999']:
    print("  q=%s  |J_quad - J_closed| = %s" % (qq, mp.nstr(abs(J_num(qq) - J_cf(qq)), 3)))
