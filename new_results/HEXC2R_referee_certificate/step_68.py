
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 40
from fractions import Fraction as Fr
import reflib as R

def qadd(x, y): return (x[0] + y[0], x[1] + y[1])
def qmul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c + b*d)
def qpow(x, n):
    out = (Fr(1), Fr(0))
    for _ in range(n): out = qmul(out, x)
    return out
ZETA = (Fr(0), Fr(1))
MU6 = [qpow(ZETA, k) for k in range(6)]
ZERO = (Fr(0), Fr(0))
def poly_mul(P, Q):
    out = [ZERO] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        for j, b in enumerate(Q):
            out[i + j] = qadd(out[i + j], qmul(a, b))
    return out
def poly_add(P, Q):
    n = max(len(P), len(Q))
    out = [ZERO] * n
    for i in range(n):
        a = P[i] if i < len(P) else ZERO
        b = Q[i] if i < len(Q) else ZERO
        out[i] = qadd(a, b)
    return out
def poly_scale(P, c):
    return [qmul(a, c) for a in P]

# replicate step_10 form_to_poly for {(4,): 1, (2,): -1}
lettervec = {(4,): Fr(1), (2,): Fr(-1)}
terms = [(c, w[0]) for w, c in lettervec.items()]
print("terms:", terms)
den = [(Fr(1), Fr(0))]
for c, j in terms:
    den = poly_mul(den, [(Fr(1), Fr(0)), (-MU6[j][0], -MU6[j][1])])
num = [ZERO]
for c, j in terms:
    p = [(Fr(1), Fr(0))]
    for c2, j2 in terms:
        if j2 == j: continue
        p = poly_mul(p, [(Fr(1), Fr(0)), (-MU6[j2][0], -MU6[j2][1])])
    cc = MU6[j]
    p = poly_scale(p, qmul((Fr(c.numerator, c.denominator), Fr(0)), cc))
    num = poly_add(num, p)
print("num:", num)
print("den:", den)
