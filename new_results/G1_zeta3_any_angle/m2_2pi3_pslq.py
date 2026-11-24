"""m2_2pi3_pslq.py -- closed form of M2 at alpha = 2pi/3 (c = 1/2) by PSLQ on the weight-3 part only.

Exact (PROVED, m2_structure.py + result.md Thm C):  M2 = z2/2 - Li2(c^2)/2 + W3,  so at c = 1/2
  W3(2pi/3) = M2 - pi^2/12 + Li2(1/4)/2.
W3 is a real weight-3 period of MT(Z[1/6]) (HPLs with letters 0, +-1 at x = 1/2).  That space has dimension 11:
 4 cubic log monomials, Li2(1/3)*{log2, log3}, pi^2*{log2, log3}, zeta(3), and two Lyndon trilogarithms.
Basis A: Lyndon pair {Li3(1/3), Li3(-1/2)};  Basis B: {Li3(2/3), Li3(-1/2)} (to show basis dependence of
the zeta(3) coefficient);  Basis C: {Li3(1/3), Li3(1/4)}.  PSLQ at 300 digits, maxcoeff 10^12.
"""
import mpmath as mp
from fractions import Fraction as Fr

mp.mp.dps = 330
M2 = None
for line in open('values_2pi_3.txt'):
    if line.startswith('M2 ='):
        M2 = mp.mpf(line.split('=')[1].strip())
W3 = M2 - mp.pi**2/12 + mp.polylog(2, mp.mpf(1)/4)/2
l2, l3, pi2 = mp.log(2), mp.log(3), mp.pi**2
common = [('zeta3', mp.zeta(3)), ('pi2*log2', pi2*l2), ('pi2*log3', pi2*l3), ('log2^3', l2**3),
          ('log2^2*log3', l2**2*l3), ('log2*log3^2', l2*l3**2), ('log3^3', l3**3),
          ('Li2(1/3)*log2', mp.polylog(2, mp.mpf(1)/3)*l2), ('Li2(1/3)*log3', mp.polylog(2, mp.mpf(1)/3)*l3)]
bases = {
    'A': [('Li3(1/3)', mp.polylog(3, mp.mpf(1)/3)), ('Li3(-1/2)', mp.polylog(3, -mp.mpf(1)/2))],
    'B': [('Li3(2/3)', mp.polylog(3, mp.mpf(2)/3)), ('Li3(-1/2)', mp.polylog(3, -mp.mpf(1)/2))],
    'C': [('Li3(1/3)', mp.polylog(3, mp.mpf(1)/3)), ('Li3(1/4)', mp.polylog(3, mp.mpf(1)/4))],
}
print('W3(2pi/3) =', mp.nstr(W3, 50))
for name, ly in bases.items():
    B = ly + common
    vec = [W3] + [v for _, v in B]
    mp.mp.dps = 300
    rel = mp.pslq(vec, maxcoeff=10**12, maxsteps=10**7)
    mp.mp.dps = 330
    if rel is None:
        print(f'basis {name}: no relation'); continue
    resid = abs(mp.fsum(r*v for r, v in zip(rel, vec)))
    terms = [(B[i][0], Fr(-rel[i + 1], rel[0])) for i in range(len(B)) if rel[i + 1] != 0]
    print(f'basis {name}: W3 = ' + ' + '.join(f'({c})*{n}' for n, c in terms))
    print(f'   integer relation: max|coef| = {max(abs(r) for r in rel)}, residual = {mp.nstr(resid, 3)}')
    z = dict(terms).get('zeta3', 0)
    print(f'   => coefficient of zeta(3) in basis {name}: {z}')
