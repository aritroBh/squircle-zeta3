"""RSTSR referee: own Arb integrators (rigorous enclosures) for M_k^main and S_k on the ray s = r*omega.
H(s) = (1/2) B(s/2,3/2) = (sqrt(pi)/4) Gamma(s/2)/Gamma(s/2+3/2)  (meromorphic, Arb-safe).
Pieces [0,r0] and [R,inf) are bounded by hand and added as error radius:
  on the ray |H| <= 1/Re s = sqrt2/r,  |E| = e^{-a r}, a = pi/(2 sqrt 2);
  M_main: |Q_k| <= r^k (r<=r0<1, since |1+s|>=1) and |Q_k|<=1/|1+s|<=1;
  S_k:    |A2| <= 2|s|^{k-1}/|s+2|^k <= 2 r^{k-1}/2^k  and  <= 2/r.
"""
from flint import arb, acb, ctx
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
def setup(prec):
    ctx.prec = prec
def consts():
    pi = arb.pi(); om = acb(arb(2).sqrt()/2, arb(2).sqrt()/2)
    return pi, om, pi/(2*arb(2).sqrt())
def H(s):
    return arb.pi().sqrt()/4*(s/2).gamma()*(s/2+arb(3)/2).rgamma()
def _int(f, r0, R, nsplit, atol):
    pi, om, a = consts()
    pts = [r0]
    while pts[-1] < 1: pts.append(pts[-1]*2)            # geometric near the poles at s=0, -1
    while pts[-1] < 24: pts.append(pts[-1] + 1)          # unit steps away from the poles
    pts += [pts[-1] + (R-pts[-1])*j/nsplit for j in range(1, nsplit+1)]
    tot = acb(0)
    for u, v in zip(pts[:-1], pts[1:]):
        tot += acb.integral(lambda r, _: f(r*om)*om, u, v, rel_tol=arb("1e-36"), abs_tol=atol, eval_limit=10**6)
    return tot
def Zmain(k, R=None, r0=None, nsplit=16, atol=None):
    """int_ray H E Q_k ds  (complex ball), Q_k = s^k/(1+s)^{k+1}."""
    pi, om, a = consts(); N = k+1
    if r0 is None: r0 = arb(2)**-7       # exact dyadic start point
    if R is None: R = arb(int((2*float((pi*N).sqrt()) + 120)/float(a)) + 10)   # e^{-aR} tiny
    f = lambda s: H(s)*(acb(0, 1)*pi*s/2).exp()*s**k/(1+s)**N
    if atol is None: atol = (-(3*pi*N).sqrt()).exp()*arb("1e-20")   # far below |U_k|
    I = _int(f, r0, R, nsplit, atol)
    err = arb(2).sqrt()*r0**k/k + 2*arb(2).sqrt()/R*(-a*R).exp()/a     # head + tail (|H|<=sqrt2/r, |Q|<=1)
    return I + acb(arb(0, err), arb(0, err))
def Z2(k, R=None, r0=None, nsplit=16, atol=None):
    """Z_k = int_ray H E A2_k ds, A2_k = -2(s+1)s^{k-1}/(s+2)^{k+1} (closed form under test)."""
    pi, om, a = consts(); N = k+1
    if r0 is None: r0 = arb(2)**-7       # exact dyadic start point
    if R is None: R = arb(int((2*float((2*pi*N).sqrt()) + 120)/float(a)) + 10)
    f = lambda s: H(s)*(acb(0, 1)*pi*s/2).exp()*(-2)*(s+1)*s**(k-1)/(s+2)**N
    if atol is None: atol = (-(3*pi*N).sqrt()).exp()*arb("1e-20")   # far below |U_k|
    I = _int(f, r0, R, nsplit, atol)
    err = 2*arb(2).sqrt()*r0**(k-1)/((k-1)*arb(2)**k) + 2*arb(2).sqrt()*2/R*(-a*R).exp()/a
    # head: int_0^r0 (sqrt2/r)*2 r^{k-1}/2^k dr = 2sqrt2 r0^{k-1}/((k-1)2^k)
    return I + acb(arb(0, err), arb(0, err))
def load_M():
    out = {}
    for line in open(os.path.join(HERE, "..", "RSTS_rest_sharpness", "mk_deg9.tsv")):
        p = line.rstrip("\n").split("\t")
        if p and p[0].isdigit():
            out[int(p[0])] = arb(p[1])
    return out
