
import os
os.chdir('research/08_theory/PRV1_depth2_proofs')
header = r'''"""
verify_full.py — PRV1: complete exact verifier for the depth-2 reductions of the
regular octagon (CX1/OCT6) and the regular decagon (DEC1).

Every step is EXACT: sympy formal differentiation algebra with exact constant pins,
or exact (finite-sum) series identities.  The high-precision numeric checks in Part 4
are independent evidence, not part of the proof.

  Part 0  Foundations: partial derivatives of Li_{2,1}; Li_{1,1} closed form.
  Part 1  Classical dilog/trilog function identities (derivative 0 + exact pin).
  Part 2  Octagon:  G, H families, derivative identities, constants, specialisation.
  Part 3  Decagon:  G_dec, H_dec, derivative identities, constants, golden lattice,
          specialisation matching DEC1 exactly.
  Part 4  Numeric evidence (200 dps) for every identity.
Run:  python verify_full.py   -> prints ALL CHECKS PASSED.
"""

from pathlib import Path
import hashlib
import pickle
import sys

import mpmath as mp
import sympy as sp

SRC_OCT = Path('research/08_theory/OCT6_octagon_closed')

def depth2(x, y, p=2, details=False):
    """Li_{p,1}(x,y) = sum_{n>m>=1} x^n y^m / (n^p m); tail <= rho^(N+1)/(1-rho)."""
    rho = max(abs(x), abs(x * y))
    if not 0 < rho < 1:
        if rho == 0:
            return (mp.mpf(0), 0, mp.mpf(0)) if details else mp.mpf(0)
        raise ValueError("The defining series requires |x|, |xy| < 1")
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target * (1 - rho)) / mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, nmax + 1):
        total += h / n**p
        h = x * h + x * xy_power / n
        xy_power *= x * y
    return (total, nmax, rho**(nmax+1)/(1-rho)) if details else total


# ======================================================================
def part0():
    print("=" * 74)
    print("PART 0: foundations (exact)")
    x, y = sp.symbols('x y')
    # -- 0.1 partial_y and partial_x as EXACT finite-sum identities (truncation N).
    N = 25
    S = sum(x**n * y**m / (n**2 * m) for n in range(1, N) for m in range(1, n))
    lhs_y = sp.diff(S, y)
    rhs_y = sp.expand((sum(x**n/n**2 for n in range(1, N))
                       - sum((x*y)**n/n**2 for n in range(1, N))/y)/(1-y))
    r = sp.expand(lhs_y - rhs_y)
    print("  0.1a d_y Li_21 truncated-series identity: residual =", r)
    assert r == 0
    lhs_x = sp.diff(S, x)
    rhs_x = sp.expand(sum(x**n * y**m / (n * m) for n in range(1, N) for m in range(1, n))/x)
    r = sp.expand(lhs_x - rhs_x)
    print("  0.1b d_x Li_21 truncated-series identity: residual =", r)
    assert r == 0
    # -- 0.2 Li_{1,1} closed form:  F = log(1-x)log(1-y) + Li2(-y/(1-y)) - Li2(-y(1-x)/(1-y))
    #     d/dx F = -log(1-x y)/(1-x) (the series derivative of Li_{1,1}); F(0,y)=0.
    z = -y*(1-x)/(1-y)
    dz = sp.diff(z, x)
    dF = (-sp.log(1-y)/(1-x)
          - (-sp.log(1-z)) * dz / z)          # d/dx[-Li2(z)] = +log(1-z) z'/z
    r = sp.factor(dF - (-sp.log(1-x*y)/(1-x)))
    print("  0.2  d/dx Li_11 closed form minus series derivative =", r)
    assert r == 0
    print("       F(0,y) = log(1)*log(1-y) + Li2(-y/(1-y)) - Li2(-y/(1-y)) = 0 (exact)")
    print("PART 0 OK")
'''
with open('verify_full.py','w') as f: f.write(header)
print("part0 written")
