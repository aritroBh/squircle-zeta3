# PSLQ of R (and of b_1''') against a basis of alternating MZVs of weight <= 5 (dimension 20).
import sys
from mpmath import mp, mpf, log, pi, zeta, polylog, pslq, nstr, fsum
fn = sys.argv[1]; D = int(sys.argv[2])
vals = dict(l.split() for l in open(fn))
mp.dps = D
l2 = log(2); z3, z5 = zeta(3), zeta(5); L4, L5 = polylog(4, mpf(1)/2), polylog(5, mpf(1)/2)
basis = [('1', 1), ('ln2', l2), ('pi^2', pi**2), ('ln2^2', l2**2), ('zeta3', z3), ('pi^2 ln2', pi**2*l2), ('ln2^3', l2**3),
 ('pi^4', pi**4), ('zeta3 ln2', z3*l2), ('pi^2 ln2^2', pi**2*l2**2), ('ln2^4', l2**4), ('Li4(1/2)', L4),
 ('zeta5', z5), ('pi^2 zeta3', pi**2*z3), ('pi^4 ln2', pi**4*l2), ('zeta3 ln2^2', z3*l2**2), ('pi^2 ln2^3', pi**2*l2**3),
 ('ln2^5', l2**5), ('Li4(1/2) ln2', L4*l2), ('Li5(1/2)', L5)]
for key in ('R', 'b1ppp'):
    x = mpf(vals[key])
    rel = pslq([x] + [b for _, b in basis], maxcoeff=10**12, maxsteps=10**7)
    print(key, 'at', D, 'digits:')
    if rel:
        print('  ', [(nm, c) for (nm, _), c in zip([('X', 0)] + basis, rel) if c])
        print('   residual', nstr(fsum(c*b for c, b in zip(rel, [x] + [b for _, b in basis])), 3))
    else: print('   none')
