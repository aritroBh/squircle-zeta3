"""RULQR referee, script 1: values of G(h) = Ehat(h)/h^3 at Chebyshev nodes on [0,H], computed with Arb
(certified numerical integration, python-flint acb.integral).  Independent of RULQ's J-table/series route:
the coefficients are later extracted by Chebyshev interpolation (rr2_extract.py).

Representation (derived independently in REFEREE.md Sec. 1): for p = 1/h > 1, q = Q(h)/h, r = q(1-h), eta = 1/q,
  L_q(p)/8 = 2^{-h} + Ehat(h),
  Ehat(h) = (h/r) int_0^oo e^{-h w/r} (1+e^{-w/r})^{-1-h} [ (1+e^{-w})^eta - 1 ] dw,
integrand <= eta e^{-w}; we integrate over [0,W] and add the tail bound to the error radius.
Usage: rr1_values.py FAM H D DPS OUT     FAM in lam2, lamhalf, p1 (q=p+1), p (q=p), gen (q=p+1/3+1/p), bdd (q=2+1/p; G=Ehat/h), sq (q=p^2; G=Ehat/h^5)
"""
import sys
from flint import acb, arb, ctx
fam, H, D, dps, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
ctx.dps = dps + 30
H = arb(H)
Qf = {'lam2': lambda h: arb(2), 'lamhalf': lambda h: arb(1) / 2, 'p1': lambda h: 1 + h, 'p': lambda h: arb(1),
      'gen': lambda h: 1 + h / 3 + h * h,           # q = p + 1/3 + 1/p: generic phase, b1 = 2/3
      'bdd': lambda h: (2 + h) * h,                 # bounded ruler q = 2 + 1/p (A2)
      'sq': lambda h: 1 / h}[fam]                   # q = p^2 (A4)
SPOW = {'bdd': 1, 'sq': 5}.get(fam, 3)             # G = Ehat / h^SPOW

def G(h):
    q = Qf(h) / h; r = q * (1 - h); eta = 1 / q
    hr, ir = h / r, 1 / r
    def f(w, a):
        A = (1 + (-w * ir).exp()).log(analytic=a)
        B = (1 + (-w).exp()).log(analytic=a)
        return (-hr * w - (1 + h) * A).exp() * (eta * B).expm1()
    W = int((dps + 20) * 2.3026) + 10
    I = acb.integral(f, 0, W, rel_tol=arb(10) ** (-(dps + 10)))
    val = hr * I.real
    tail = hr * eta * arb(-W).exp()            # int_W^oo eta e^{-w} dw
    val = val + arb(0, tail.upper())           # widen by the tail
    return val / h ** SPOW

pi = arb.pi()
with open(out, 'w') as fo:
    for j in range(D):
        h = H / 2 * (1 + (pi * (j + arb(1) / 2) / D).cos())
        g = G(h)
        fo.write('%s %s\n' % (h.str(dps + 20, radius=False), g.str(dps + 5, radius=True)))
        fo.flush()
