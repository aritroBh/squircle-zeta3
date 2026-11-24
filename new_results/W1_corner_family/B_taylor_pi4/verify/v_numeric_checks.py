"""v_numeric_checks.py -- verifier's numerical checks at the square corner c0 = 1/sqrt2 (own HPL code).
(a) exact: the weight-2 HPL part of G1's R_3 and the weight-0 part of Psi (=> W_2(c) = z2/2 - Li2(c^2)/2 exactly);
(b) NUMERICAL: B3(c0) vs the squircle M_3 (paper's stored 1200-digit value), and the weight split
    W_4(c0), W_3'(c0), W_2(c0) vs the note's W_3', T  (W_3' typed by the verifier from main.tex Theorem 3,
    as transcribed in the solver's check_b3_split.py -- re-read below against main.tex);
(c) the prefactor components at c0.
"""
import json, re, sys
from fractions import Fraction as Fr
import mpmath as mp
sys.path.insert(0, '.')
from v_nested_proof import load

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 160
mp.mp.dps = DPS + 25
PAPER = '<home>)/paper/squircle_note/outputs/'
R3 = load('hpl_fit_k3_(0,0,0).json'); PSI = load('hpl_fit_k3_(1,0,0).json'); MPHI = load('hpl_fit_k3_(0,1,0).json')

print('(a) weight-2 HPL terms of R_3:', [(e, w, str(v)) for e, w, v in R3 if len(w) == 2])
print('    weight-0 terms of Psi:', [(e, w, str(v)) for e, w, v in PSI if len(w) == 0],
      '; weight-1 terms of R_3:', [(e, w, str(v)) for e, w, v in R3 if len(w) == 1])
print('    weight-0/1 terms of -Phi:', [(e, w, str(v)) for e, w, v in MPHI])

c0 = 1/mp.sqrt(2)
NT = int((DPS + 30)*mp.log(10)/mp.log(mp.sqrt(2))) + 20
cache = {}
def coeffs(w):
    if w in cache:
        return cache[w]
    if len(w) == 1:
        out = [mp.mpf(0)] + [(mp.mpf(1)/m if w[0] == 1 else mp.mpf((-1)**(m + 1))/m) for m in range(1, NT + 1)]
    else:
        inner = coeffs(w[1:]); a = w[0]
        out = [mp.mpf(0)]*(NT + 1)
        if a == 0:
            for m in range(1, NT + 1):
                out[m] = inner[m]/m
        else:
            run = mp.mpf(0)
            for m in range(0, NT):
                run = run*a + inner[m]
                out[m + 1] = run/(m + 1)
    cache[w] = out
    return out
pw = [mp.mpf(1)]
for _ in range(NT):
    pw.append(pw[-1]*c0)
def H(w):
    cf = coeffs(w)
    return mp.fsum(cf[m]*pw[m] for m in range(1, NT + 1))

z2, z3 = mp.zeta(2), mp.zeta(3)
def ev(terms, wt=None, e_only=None, with_pref=True):
    s = mp.mpf(0)
    for e, w, v in terms:
        if wt is not None and len(w) != wt:
            continue
        if e_only is not None and e != e_only:
            continue
        hv = mp.mpf(1) if not w else H(w)
        s += mp.mpf(v.numerator)/v.denominator*(c0**e if with_pref else 1)*hv
    return s
M3c = ev(R3) + z2*ev(PSI) + z3*ev(MPHI)
M3s = mp.mpf(open(PAPER + 'M3_methodA_1200.txt').read().strip())
M3sB = mp.mpf(open(PAPER + 'M3_methodB_1200.txt').read().strip()) if __import__('os').path.exists(PAPER + 'M3_methodB_1200.txt') else None
print(f'(b) B3(c0) - M3_square(methodA) = {mp.nstr(M3c - M3s, 5)}' + (f';  - methodB = {mp.nstr(M3c - M3sB, 5)}' if M3sB is not None else ''))
W4 = ev(R3, 4) + z2*ev(PSI, 2) + z3*ev(MPHI, 1)
W3 = ev(R3, 3) + z2*ev(PSI, 1) + z3*ev(MPHI, 0)
W2 = ev(R3, 2) + z2*ev(PSI, 0)
r2 = mp.sqrt(2); b = r2 - 1; lam = mp.log(2); L = mp.log(1 + r2)
W3p = ((7*r2/32 - mp.mpf(3)/16)*z3 + (1 + r2)*mp.polylog(3, b) + (1 - r2)*mp.polylog(3, -b)
       - 2*r2*mp.polylog(3, 1/r2) + (2*L - r2*lam)*mp.polylog(2, b) + (r2/12 - mp.mpf(1)/8)*mp.pi**2*L
       - (1 + r2)/24*mp.pi**2*lam + (1 + r2)/24*lam**3 + r2/4*lam**2*L - (1 + r2)/2*lam*L**2 + (5 - r2)/6*L**3)
T = (mp.pi**2 + 6*lam**2)/24
print(f'    W2(c0) - T = {mp.nstr(W2 - T, 5)};  z2/2 - Li2(1/2)/2 - T = {mp.nstr(z2/2 - mp.polylog(2, 0.5)/2 - T, 5)}')
print(f"    W3'(c0) - W3'_note = {mp.nstr(W3 - W3p, 5)}")
print(f'    W4(c0) - (M3_square - W3_note - T) = {mp.nstr(W4 - (M3s - W3p - T), 5)}')
print('    W4(c0) =', mp.nstr(W4, 60))
for e in (-1, 0, 1):
    Ae = ev(R3, 4, e, False) + z2*ev(PSI, 2, e, False) + z3*ev(MPHI, 1, e, False)
    Ae_pure = ev(R3, 4, e, False)
    print(f'(c) A_{e}(c0) incl. zeta terms = {mp.nstr(Ae, 30)};  pure weight-4 HPL part only = {mp.nstr(Ae_pure, 30)}')
