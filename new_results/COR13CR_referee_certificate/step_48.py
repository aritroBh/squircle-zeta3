import time, sys
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import reflib13 as R13
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

qadd, qsub, qmul, qdiv, qneg = R13.qadd, R13.qsub, R13.qmul, R13.qdiv, R13.qneg
MU6 = R13.MU6
ZERO = (Fr(0), Fr(0)); ONE = (Fr(1), Fr(0))

def poly_mul(P, Q):
    out = [ZERO] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        if a == ZERO: continue
        for j, b in enumerate(Q):
            if b != ZERO:
                out[i+j] = qadd(out[i+j], qmul(a, b))
    return out
def poly_add(P, Q):
    n = max(len(P), len(Q))
    out = [ZERO] * n
    for i in range(n):
        a = P[i] if i < len(P) else ZERO
        b = Q[i] if i < len(Q) else ZERO
        out[i] = qadd(a, b)
    return out
def poly_scale(P, c): return [qmul(a, c) for a in P]
def trim(P):
    P = list(P)
    while len(P) > 1 and P[-1] == ZERO: P.pop()
    return P

def letterform_to_poly(lettervec):
    """sum_j c_j letter_j(u) as (num, den): letter -1 = 1/u; letter j>=0 = zeta^j/(1-zeta^j u)."""
    terms = [(c, w[0]) for w, c in lettervec.items()]
    den = [ONE]
    if any(j == -1 for c, j in terms):
        den = poly_mul(den, [ZERO, ONE])
    for c, j in terms:
        if j >= 0:
            den = poly_mul(den, [ONE, qneg(MU6[j])])
    num = [ZERO]
    for c, j in terms:
        p = [ONE]
        for c2, j2 in terms:
            if j2 == j: continue
            p = poly_mul(p, [ZERO, ONE] if j2 == -1 else [ONE, qneg(MU6[j2])])
        cc = ONE if j == -1 else MU6[j]
        p = poly_scale(p, qmul((Fr(c), Fr(0)), cc))
        num = poly_add(num, p)
    return trim(num), trim(den)

def dlog_linear(alpha, beta):
    """dlog(alpha u + beta) = alpha/(alpha u + beta) du -> (num, den)."""
    return [alpha], trim([beta, alpha])

def rat_sub(N1, D1, N2, D2):
    """N1/D1 - N2/D2 -> (num, den)."""
    num = poly_add(poly_mul(N1, D2), poly_scale(poly_mul(N2, D1), (Fr(-1), Fr(0))))
    den = poly_mul(D1, D2)
    return trim(num), trim(den)

def rat_eq(N1, D1, N2, D2):
    """cross-multiply comparison of rational functions over Q(zeta)."""
    return trim(poly_mul(N1, D2)) == trim(poly_mul(N2, D1))

ok = 0; total = 0
for mi, E in enumerate(R13.maps13()):
    f0, fh, f1, fm, finf = [MU6[e] for e in E]
    # f: A=finf, C=1, D=(finf-f1)/(f1-f0), B=f0*D ; g=f^{-1}: A'=D, B'=-B, C'=-1, D'=A=finf
    D = qdiv(qsub(finf, f1), qsub(f1, f0))
    B = qmul(f0, D)
    Ap, Bp, Cp, Dp = D, qneg(B), qneg(ONE), finf
    for sig, q in ((0, (Fr(0), Fr(0))), (1, ONE), (-1, qneg(ONE))):
        total += 1
        # LHS: dlog(g(u)-q) = dlog((Ap-qCp)u + (Bp-qDp)) - dlog(Cp u + Dp)
        A2 = qsub(Ap, qmul(q, Cp)); B2 = qsub(Bp, qmul(q, Dp))
        N1, D1 = dlog_linear(A2, B2)
        N2, D2 = dlog_linear(Cp, Dp)
        NL, DL = rat_sub(N1, D1, N2, D2)
        # RHS: my pullback vector; for sig=1 the HPL form is -dlog(x-1)
        vec = R13.hpl_form_pullback13(mi, sig)
        NR, DR = letterform_to_poly(vec)
        assert rat_eq(NL, DL, NR, DR), (mi, sig)
        ok += 1
log(f"exact 1-form pullback identities verified over Q(zeta): {ok}/{total} (12 maps x 3 HPL letters)")
assert ok == total
with open('progress.log','a') as f:
    f.write("referee step(forms exact): 36 pulled-back 1-form identities verified EXACTLY over Q(zeta) (12 maps x 3 letters)\n")
log("STEP PASS")
