# 13_proof_Q0.py -- exact evaluation of the vertex constant Q0 = 6 (P1 - P2a - P2c).
# Each piece is reduced, by explicit antiderivatives, to real dilogarithms at S-units of Q(sqrt2)
# plus elementary terms; complex dilogarithms are removed by reflection / Landen / inversion and
# Lemma A:  Re Li2(b e^{+-i pi/4}) = Li2(b)/2 + pi^2/192.  The remaining real Li2 values are
# reduced to the basis {Li2(b), Li2(-b^4), Li2(1/3)} by dilog_engine.py (exact linear algebra on
# instances of classical functional equations).  Every intermediate step is also checked numerically.
import mpmath as mp, json
from fractions import Fraction as Fr
from qs2 import Q2, W2, Li2, PI2, PIA, R2, A, B, logprod, logsq, logvec
from dilog_engine import pool, Reducer
mp.mp.dps = 60
ONE = Q2(1)
def chk(name, w, val):
    d = w.num() - val
    print("  %-38s numeric check %s" % (name, mp.nstr(d, 3)))
    assert abs(d) < mp.mpf(10)**-45, name

# ---------- generic real primitive  int_A^B log(v+p) dv/(v+r),  v+p>0 on [A,B], v+r of constant sign
def I(p, r, Aa, Bb):
    p, r, Aa, Bb = Q2(1)*p, Q2(1)*r, Q2(1)*Aa, Q2(1)*Bb
    w = W2()
    if p == r:
        return logsq(Bb+p, Fr(1, 2)) + logsq(Aa+p, Fr(-1, 2))
    c = p - r
    ylo, yhi = Aa + r, Bb + r
    if float(ylo) > 0:
        if float(c) > 0:   # log(y+c) = log c + log(1+y/c)
            return (logprod(c, yhi) - logprod(c, ylo) - Li2(-yhi/c) + Li2(-ylo/c))
        d = -c             # y > d > 0: log(y-d) = log y + log(1-d/y)
        return logsq(yhi, Fr(1, 2)) - logsq(ylo, Fr(1, 2)) + Li2(d/yhi) - Li2(d/ylo)
    else:                  # y<0: Y=-y>0, log(c-Y) = log c + log(1-Y/c), dy/y = dY/Y
        Yhi, Ylo = -yhi, -ylo
        return (logprod(c, Yhi) - logprod(c, Ylo) - Li2(Yhi/c) + Li2(Ylo/c))
def Iconst(cst, r, Aa, Bb):   # log(cst) * int dv/(v+r)
    r = Q2(1)*r
    return logprod(cst, (Q2(1)*Bb + r)) - logprod(cst, (Q2(1)*Aa + r))

def num_int(f, a_, b_): return mp.quad(f, [a_, b_])
r2 = mp.sqrt(2); a = 1 + r2; b = r2 - 1; lam = mp.log(2); L = mp.log(a); l3 = mp.log(3)
alpha = mp.atan(1/r2)

# Lemma A (proved by differentiation in theta; checked here)
ReLi2_bz = Li2(B, Fr(1, 2)) + PI2(Fr(1, 192))
chk("Lemma A: Re Li2(b e^{i pi/4})", ReLi2_bz, mp.re(mp.polylog(2, b*mp.expjpi(mp.mpf(1)/4))))
# general family check at a few theta
for th in [mp.mpf('0.3'), mp.mpf('0.9')]:
    rr = 2*mp.cos(th) - 1
    lhs = mp.re(mp.polylog(2, rr*mp.exp(1j*th))); rhs = mp.polylog(2, rr)/2 + 3*th**2/4 - mp.pi*th/2 + mp.pi**2/12
    print("  Lemma A family at theta=%s: %s" % (th, mp.nstr(lhs - rhs, 3)))

# ================================================================ P1
qv = lambda v: (v*v+1)/(2*r2*v)
print("P1:")
# complex part: 2C(a) - 2C(b) + sqrt2 C(0)
z1 = (1 + 1j*b)/2; z2 = 1 - a*mp.expjpi(mp.mpf(1)/4); z3 = -(b + 1j)/2
ReLi2_z1 = PI2(Fr(1, 12)) + logsq(B/R2, Fr(-1, 8)) + PI2(Fr(-1, 128))          # Re z=1/2: pi^2/12 - (log^2|z|+arg^2)/2
chk("Re Li2(z1)", ReLi2_z1, mp.re(mp.polylog(2, z1)))
ReLi2_z2 = (PI2(Fr(1, 3)) + ReLi2_bz + logsq(A, Fr(1, 2)) + PI2(Fr(-9, 32))
            - logprod(R2*A, A, Fr(1, 2)) + PI2(Fr(-5, 32)))
chk("Re Li2(z2)", ReLi2_z2, mp.re(mp.polylog(2, z2)))
ReLi2_z3 = -ReLi2_bz + logsq(A/R2, Fr(-1, 8)) + PI2(Fr(1, 128))
chk("Re Li2(z3)", ReLi2_z3, mp.re(mp.polylog(2, z3)))
Ca = ReLi2_z1.scale(2)
Cb = ReLi2_z2.scale(-2) + ReLi2_z3.scale(2)
C0 = Li2(-(A*A), Fr(-1, 2)) + PI2(Fr(-1, 24))
chk("C(a)", Ca, num_int(lambda v: (mp.log(v*v+1)-mp.log(a*a+1))/(v-a), 1, a))
chk("C(b)", Cb, num_int(lambda v: (mp.log(v*v+1)-mp.log(b*b+1))/(v-b), 1, a))
chk("C(0)", C0, num_int(lambda v: mp.log(v*v+1)/v, 1, a))
# real parts
Ta_real = Li2(B/R2, -2) + Li2(Q2(Fr(1, 2)), -2)
chk("T_a real part", Ta_real, 2*num_int(lambda v: (mp.log(4*r2*a)-mp.log(v+a)-mp.log(v+b))/(v-a), 1, a))
Tb_real = (I(B, -B, 1, A) - Iconst(2*B, -B, 1, A) + I(A, -B, 1, A) - Iconst(A+B, -B, 1, A)).scale(2)
chk("T_b real part", Tb_real, -2*num_int(lambda v: (mp.log(4*r2*b)-mp.log(v+a)-mp.log(v+b))/(v-b), 1, a))
Tc = W2()
for p, cp in ((A, 1), (B, 1), (0, -1)):
    for rr, cr in ((B, 1), (A, -1)):
        Tc = Tc + (I(p, rr, 1, A) if p != 0 else I(Q2(0), rr, 1, A)).scale(-2*cp*cr)
for rr, cr in ((B, 1), (A, -1)):
    Tc = Tc + Iconst(2*R2, rr, 1, A).scale(2*cr)
chk("T_c", Tc, -2*num_int(lambda v: mp.log(1+qv(v))*(1/(v+b)-1/(v+a)), 1, a))
Td_real = (logprod(R2, A) + logsq(A, Fr(1, 2))).scale(-R2)
P1 = Ca.scale(2) - Cb.scale(2) + C0.scale(R2) + Ta_real + Tb_real + Tc + Td_real
d = json.load(open("pieces_300.json"))
chk("P1 (exact, unreduced)", P1, mp.mpf(d["P1"]))

# ================================================================ P2c
print("P2c:")
ReChi2_bz = ReLi2_bz - Li2(-(B**4), Fr(1, 16))
chk("Re chi2(b e^{-i pi/4})", ReChi2_bz, mp.re((mp.polylog(2, b*mp.expjpi(-mp.mpf(1)/4)) - mp.polylog(2, -b*mp.expjpi(-mp.mpf(1)/4)))/2))
x = (r2 + 1j)/3
ReChi2_x = PI2(Fr(1, 8)) + logprod(Q2(3), A, Fr(-1, 4)) + PIA(Fr(-1, 8)) - ReChi2_bz
chk("Re chi2((sqrt2+i)/3)  [Landen]", ReChi2_x, mp.re((mp.polylog(2, x) - mp.polylog(2, -x))/2))
P2c_over = ReChi2_x.scale(-2) + PI2(Fr(1, 4)) + PIA(Fr(-1, 2))
P2c = P2c_over.scale(R2)
chk("P2c (exact)", P2c, mp.mpf(d["P2c"]))

# ================================================================ P2a
print("P2a:")
# W = V/4 + (Li2(2/3)-Li2(1/2))/4,  V = int_0^{1/sqrt2} [2log(1+t)-log(sqrt2-t)-log(sqrt2+t)][1/(sqrt2-t)-1/(sqrt2+t)] dt
lo, hi = Q2(0), Q2(0, Fr(1, 2))     # 1/sqrt2 = sqrt2/2
V = W2()
for p, cp in ((ONE, 2), (-R2, -1), (R2, -1)):       # log(t+p): t+1, sqrt2-t = -(t-sqrt2) handled below
    for rr, cr in ((-R2, -1), (R2, -1)):            # 1/(sqrt2-t) = -1/(t-sqrt2);  -1/(sqrt2+t)
        if p == -R2:
            # log(sqrt2 - t): substitute s=-t (t in [0,1/sqrt2] -> s in [-1/sqrt2,0]); int log(s+sqrt2) ds/(s - rr) ... use
            # int_0^{h} log(sqrt2-t) dt/(t+rr) = int_{-h}^{0} log(sqrt2+s) ds/(s - rr)
            term = I(R2, -rr, -hi, lo).scale(-1)
        else:
            term = I(p, rr, lo, hi)
        V = V + term.scale(cp*cr)
Vnum = num_int(lambda t: (2*mp.log(1+t)-mp.log(r2-t)-mp.log(r2+t))*(1/(r2-t)-1/(r2+t)), 0, 1/r2)
chk("V", V, Vnum)
Wexact = V.scale(Fr(1, 4)) + Li2(Q2(Fr(2, 3)), Fr(1, 4)) + Li2(Q2(Fr(1, 2)), Fr(-1, 4))
phi0 = mp.asinh(1/r2)
chk("W = int arccosh(sqrt(2/3)cosh phi)", Wexact, mp.quad(lambda p: mp.acosh(mp.sqrt(mp.mpf(2)/3)*mp.cosh(p)), [phi0, L]))
P2a = (logprod(Q2(3), A, Fr(1, 2)) - Wexact).scale(R2)
chk("P2a (exact)", P2a, mp.mpf(d["P2a"]))

# ================================================================ reduction
open("outputs_13_certificates.txt", "w").close()
print("reduction to the basis {Li2(b), Li2(-b^4), Li2(1/3)}:")
P = pool()
red = Reducer(P, basis=[B, -(B**4), Q2(Fr(1, 3))])
out = {}
for name, w in (("P1", P1), ("P2a", P2a), ("P2c", P2c)):
    missing = [k for k in w.li2keys() if k not in P]
    if missing: print("  WARNING points not in pool:", missing)
    r, certs = red.reduce(w)
    out[name] = r
    with open("outputs_13_certificates.txt", "a") as fh:
        fh.write("=== %s: functional-equation instances used (multiplier * instance), rational part / sqrt2 part\n" % name)
        for part, cert in zip(("1", "sqrt2"), certs):
            for k, v in sorted(cert.items(), key=lambda kv: str(kv[0])):
                fh.write("  [%s] %s * %s\n" % (part, v, k))
    print("  %s = %r" % (name, r))
    chk("  reduced %s" % name, r, mp.mpf(d[name]))
Q0 = (out["P1"] - out["P2a"] - out["P2c"]).scale(6)
Q0r, _ = red.reduce(Q0)
print("Q0 =", Q0r)
chk("Q0", Q0r, mp.mpf(d["Q0"]))
# target closed form
target = (PI2(Fr(-1, 8)) + PI2(Q2(0, Fr(1, 16))) + PIA(Q2(0, Fr(3, 2))) + logsq(Q2(2), -3)
          + logprod(Q2(3), A, Q2(0, -9)) + logsq(A, 6) + logsq(A, Q2(0, 6)) + Li2(B, Q2(0, -6))
          + Li2(-(B**4), Q2(0, Fr(9, 4))))
diff, _ = red.reduce(Q0r - target)
print("Q0 - [closed form of Theorem]  (exact) =", diff if diff.d else "0")
