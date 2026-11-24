# g9_hex5_compare.py -- Theorem 26.  (1) compares the table of 307 monomials rebuilt by g9_hex5_chain.sh (directory $G9WORK)
# with the table of the notes (K1g_hexagon_order5_assembly/coef_acc_raw.pkl): exact equality of every coefficient;
# (2) census of the table; (3) the coefficients of zeta(5) and zeta(2)zeta(3) against Z_e^(5), Z_e^(2,3) of Proposition 9;
# (4) writes the plain-text table g9_hex5_table.txt (E_e = 120 * sum of the lines) and prints its sha256.
import os, sys, pickle, hashlib, time, sympy as sp
WD = os.environ['G9WORK']
NOTES = '<project>/research/08_theory/K1g_hexagon_order5_assembly/coef_acc_raw.pkl'
t0 = time.time()
mine = pickle.load(open(f'{WD}/coef_acc_raw.pkl', 'rb'))
theirs = pickle.load(open(NOTES, 'rb'))
print('monomials: rebuilt', len(mine), ' notes', len(theirs), ' same key set:', set(mine) == set(theirs))
bad = 0
for k in mine:
    if sp.expand(mine[k] - theirs[k]) != 0: bad += 1
print('coefficients differing (exact, sympy expand of the difference):', bad)
# census
from collections import Counter
xp = Counter(k[1] for k in mine); print('x-exponents:', dict(sorted(xp.items())))
atoms = Counter(a for k in mine for (a, e) in k[0]); print('atoms occurring:', len(atoms), sorted(atoms))
zs = Counter(k[2] for k in mine); print('(z2,z3,z5)-exponents occurring:', dict(sorted(zs.items())))
print('monomials with z5 or with z2*z3 (z2>=1 and z3>=1):', [(k[0], k[1], k[2]) for k in mine if k[2][2] >= 1 or (k[2][0] >= 1 and k[2][1] >= 1)])
# powers of sqrt(L2)
L2 = sp.Symbol('L2')
pw = Counter()
for k, c in mine.items():
    for term in sp.Add.make_args(sp.expand(c)):
        pw[term.as_powers_dict().get(L2, 0)] += 1
print('powers of L2 in the coefficients (L2^e, e half-integers):', dict(sorted(pw.items())), ' -> sqrt(L2) = L, so L^-9 ... L^1')
# weight: (sum of word lengths of the atoms) + 2 a + 3 b + 5 c + (degree in lA0, lA1, lB) <= 5 for every term of every coefficient
WLEN = {'H1': 1, 'Li2': 2, 'Li3': 3, 'H011': 3, 'Li4': 4, 'H0011': 4, 'H0111': 4, 'Li5': 5, 'H00011': 5, 'H00101': 5, 'H00111': 5, 'H01011': 5, 'H01111': 5}
lA0s, lA1s, lBs = sp.symbols('lA0 lA1 lB')
maxw = 0
for k, c in mine.items():
    w0 = sum(WLEN[a]*e for a, e in k[0]) + 2*k[2][0] + 3*k[2][1] + 5*k[2][2]
    for term in sp.Add.make_args(sp.expand(c)):
        pdict = term.as_powers_dict()
        maxw = max(maxw, w0 + sum(pdict.get(l, 0) for l in (lA0s, lA1s, lBs)))
print('largest weight of a term of the table (atoms + zeta weights + degree in lA0,lA1,lB):', maxw)
# zeta(5), zeta(2)zeta(3)
A0, A1, B, kap, q0, V02 = sp.symbols('A0 A1 B kap q0 V02')
rel = {A1: 1 - B, A0: 1 - B - B*kap}
L = sp.sqrt(L2)
xi0, xi1 = A0, A1
L0 = -q0/L; L1 = L + q0/L; h2 = V02 - q0**2/L2
xb = (xi0*L0 + xi1*L1)/L; xp_, xt = xi0 + xi1, xi0*xi1; LL = L0*L1/L2; m = h2/L2
Z5 = 360*xt*(h2/L)*(20*xb**3 - 32*xb**2 + 15*xb - 2 - 4*(5*xb - 4)*xt*LL + m*(5*xt*xp_ + 8*xp_**2 - 28*xt - 15*(xp_**2 - 3*xt)*xb))
Z23 = 120*xt*(h2/L)*(12*xb**2 - 12*xb + 2 - 12*(5*xb - 4)*xt*LL + m*(15*xt*xp_ - 3*xp_**2 + 24*xt - 45*xt*xb))
k5 = ((), 0, (0, 0, 1)); k23 = ((), 0, (1, 1, 0))
d5 = sp.simplify((120*mine[k5] - Z5).subs(rel)); d23 = sp.simplify((120*mine[k23] - Z23).subs(rel))
print('120 x [zeta5] - Z_e^(5) =', d5, ';  120 x [zeta2 zeta3] - Z_e^(2,3) =', d23, ' (after A1 = 1-B, A0 = 1-B-B kappa)')
# plain-text table
out = []
for k in sorted(mine, key=lambda k: (k[1], k[2], k[0])):
    atoms_s = '*'.join('%s^%d' % (a, e) for a, e in k[0]) or '1'
    out.append('x^%d | %s | z2^%d z3^%d z5^%d | %s' % (k[1], atoms_s, k[2][0], k[2][1], k[2][2], str(sp.factor(mine[k])) if False else str(sp.expand(mine[k]))))
txt = '# Theorem 26. Table of the 307 monomials.  E_e^(5) = 120 * sum over lines of  coefficient * x^k * (product of the atoms) * z2^a z3^b z5^c,\n' \
      '# with z_i = zeta(i), x = 1 - r = -kap, atoms = H1 (= l), Li2, Li3, H011, Li4, H0011, H0111, Li5, H00011, H00101, H00111, H01011, H01111 (harmonic polylogarithms at x),\n' \
      '# coefficients polynomial in A0, A1, B, kap, q0, V02, lA0, lA1, lB and in L2 = L^2 (with sqrt(L2) = L), where A1 = 1 - B, A0 = 1 - B - B kap.\n' \
      '# format:  x-power | atoms | zeta powers | coefficient\n' + '\n'.join(out) + '\n'
here = os.path.dirname(os.path.abspath(__file__))
open(os.path.join(here, 'g9_hex5_table.txt'), 'w').write(txt)
print('wrote g9_hex5_table.txt: %d lines, %d bytes, sha256 %s' % (len(out), len(txt.encode()), hashlib.sha256(txt.encode()).hexdigest()))
print('done (%.0f s)' % (time.time()-t0))
