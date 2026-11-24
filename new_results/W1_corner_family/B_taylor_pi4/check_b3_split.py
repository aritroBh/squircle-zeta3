"""check_b3_split.py -- NUMERICAL.
Evaluate G1's HPL form B3 of M_3(c) (data read from G1's hpl_fit_k3_*.json, read-only) at c = 1/sqrt2,
compare with the squircle M_3 (1200-digit stored value), and split B3 by HPL weight:
   M_3(c) = W_4(c) + W_3'(c) + W_2(c),
   W_4(c) = [weight-4 HPL terms of R_3] + z2 * [weight-2 part of Psi] + z3 * [weight-1 part of -Phi],
then compare W_4(1/sqrt2) with the squircle W_4 = M_3 - W_3' - T (note's Theorem 3).
"""
import json, re, sys
from fractions import Fraction as Fr
import mpmath as mp
sys.path.insert(0, '.')
import hpl

G1 = 'research/10_frontier/G1_zeta3_any_angle/'
PAPER = '<home>)/paper/squircle_note/outputs/'


def load(name):
    out = []
    for nm, v in json.load(open(G1 + name)):
        if nm == '1':
            out.append((0, (), Fr(v)))
            continue
        if nm == 'c':
            out.append((1, (), Fr(v)))
            continue
        m = re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(2).split(','))
        out.append((e, w, Fr(v)))
    return out


R3 = load('hpl_fit_k3_(0,0,0).json')
PSI = load('hpl_fit_k3_(1,0,0).json')
MPHI = load('hpl_fit_k3_(0,1,0).json')   # coefficient of z3 = -Phi


def ev(terms, c, N, weight=None):
    s = mp.mpf(0)
    for e, w, v in terms:
        if weight is not None and len(w) != weight:
            continue
        hv = mp.mpf(1) if w == () else hpl.value(w, c, N)
        s += mp.mpf(v.numerator)/v.denominator * c**e * hv
    return s


if __name__ == '__main__':
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 110
    mp.mp.dps = dps + 20
    c = 1/mp.sqrt(2)
    N = int((dps + 25)*mp.log(10)/mp.log(mp.sqrt(2))) + 10
    z2, z3 = mp.zeta(2), mp.zeta(3)
    M3 = ev(R3, c, N) + z2*ev(PSI, c, N) + z3*ev(MPHI, c, N)
    M3s = mp.mpf(open(PAPER + 'M3_methodA_1200.txt').read().strip())
    print('B3 at c=1/sqrt2 minus squircle M3 (methodA):', mp.nstr(M3 - M3s, 5))
    W4c = ev(R3, c, N, 4) + z2*ev(PSI, c, N, 2) + z3*ev(MPHI, c, N, 1)
    W3c = ev(R3, c, N, 3) + z2*ev(PSI, c, N, 1) + z3*ev(MPHI, c, N, 0)
    W2c = ev(R3, c, N, 2) + z2*ev(PSI, c, N, 0)
    print('weight split check (sum - M3):', mp.nstr(W4c + W3c + W2c - M3, 5))
    # squircle pieces (note, Theorem 3 and Prop. T)
    r2 = mp.sqrt(2); b = r2 - 1; lam = mp.log(2); L = mp.log(1 + r2)
    W3p = ((7*r2/32 - mp.mpf(3)/16)*z3 + (1 + r2)*mp.polylog(3, b) + (1 - r2)*mp.polylog(3, -b)
           - 2*r2*mp.polylog(3, 1/r2) + (2*L - r2*lam)*mp.polylog(2, b) + (r2/12 - mp.mpf(1)/8)*mp.pi**2*L
           - (1 + r2)/24*mp.pi**2*lam + (1 + r2)/24*lam**3 + r2/4*lam**2*L - (1 + r2)/2*lam*L**2
           + (5 - r2)/6*L**3)
    T = (mp.pi**2 + 6*lam**2)/24
    W4s = M3s - W3p - T
    print('W2(1/sqrt2) - T          :', mp.nstr(W2c - T, 5))
    print("W3'(1/sqrt2) - W3'_square:", mp.nstr(W3c - W3p, 5))
    print('W4(1/sqrt2) - W4_square  :', mp.nstr(W4c - W4s, 5))
    print('W4(1/sqrt2) =', mp.nstr(W4c, 60))
    # pieces by prefactor c^e of the weight-4 part
    for e in (-1, 0, 1):
        part = sum((mp.mpf(v.numerator)/v.denominator*hpl.value(w, c, N) for (ee, w, v) in R3 if ee == e and len(w) == 4), mp.mpf(0))
        part += z2*sum((mp.mpf(v.numerator)/v.denominator*hpl.value(w, c, N) for (ee, w, v) in PSI if ee == e and len(w) == 2), mp.mpf(0))
        part += z3*sum((mp.mpf(v.numerator)/v.denominator*hpl.value(w, c, N) for (ee, w, v) in MPHI if ee == e and len(w) == 1), mp.mpf(0))
        print(f'  pure component W4^({e})(1/sqrt2) =', mp.nstr(part, 40))
