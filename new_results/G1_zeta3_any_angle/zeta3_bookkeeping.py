"""zeta3_bookkeeping.py -- how the coefficient of zeta(3) in M2 depends on the chosen trilogarithm basis.

Uniform classical formula (m2_classical_assemble.py):  M2 = z2 Phi(c) + sum_i C_i(c) Li3(a_i(c)) + ... + Z(c) zeta(3),
Z(c) = 7/4 + 2c, with a_i in {c, -c, 1-c, 1/(1+c), (1+c)/2, (1-c)/(1+c), -(1-c)/(1+c), 2c/(1+c)}.
 (1) FLAT NORMALISATION: replace each Li3(a_i(c)) by Li3(a_i(c)) - Li3(a_i(0)) (functions vanishing at the
     flat angle c=0).  The new zeta(3)-coefficient is Z(c) + sum_i C_i(c) z(a_i(0)), where z(v) is the zeta(3)
     part of Li3(v) for v in {0, 1, 1/2, -1}.  Also tracks pi^2 log2 and log^3 2.  Expected: all vanish
     identically (the canonical statement of result.md Thm A).
 (2) SQUIRCLE c = 1/sqrt2: reduce the eight Li3 values to the note's basis {Li3(b), Li3(-b), Li3(1/sqrt2), zeta3}
     using the zeta(3)-parts of the note's Lemma (lem:trilog; each identity re-checked numerically here, 60 digits)
     and D3: Li3(b^2) = 4Li3(b)+4Li3(-b).  Expected: 189 sqrt2/64 - 19/16 (note's Theorem M2).
 (3) 2pi/3, c = 1/2: reduce Li3(1/2), Li3(2/3), Li3(3/4), Li3(-1/3) to basis A = {Li3(1/3), Li3(-1/2)} by PSLQ
     (12-element bases, 200 digits).  Expected: 3/16 (m2_2pi3_pslq.py, basis A).
"""
import json, re
import sympy as sp
import mpmath as mp

c = sp.symbols('c', positive=True)
R2 = json.load(open('hpl_fit_k2_(0,0,0).json'))
CD = json.load(open('classical_dict.json'))
coef = {}
for name, cstr in R2:
    mm = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    e = int(mm.group(1)); w = [int(z) for z in mm.group(2).split(',')]
    pref = sp.Rational(cstr)*c**e
    if len(w) == 2:
        continue
    for fn, cn, fr in CD[str(w)]:
        coef[(fn, cn)] = sp.simplify(coef.get((fn, cn), 0) + pref*sp.Rational(fr))
Z = sp.factor(coef[('1', 'z3')])
C = {fn[4:-1]: sp.factor(v) for (fn, cn), v in coef.items() if fn.startswith('Li3(') and cn == '1'}
print('Z(c) =', Z)
print('Li3 coefficients C_i(c):', C)

# (1) flat normalisation: a_i(0) and zeta(3), pi^2 log2, log^3 2 content of Li3(a_i(0))
a0 = {'x': 0, '-x': 0, '1-x': 1, '1/(1+x)': 1, '(1+x)/2': sp.Rational(1, 2), '(1-x)/(1+x)': 1,
      '-(1-x)/(1+x)': -1, '2x/(1+x)': 0}
content = {0: (0, 0, 0), 1: (1, 0, 0), -1: (sp.Rational(-3, 4), 0, 0),
           sp.Rational(1, 2): (sp.Rational(7, 8), sp.Rational(-1, 12), sp.Rational(1, 6))}   # (z3, pi2 log2, log2^3)
tot = [Z, sp.factor(coef.get(('1', 'pi2*log2'), 0)), sp.factor(coef.get(('1', 'log2^3'), 0))]
for arg, Ci in C.items():
    zc = content[a0[arg]]
    for k in range(3):
        tot[k] += Ci*zc[k]
print('(1) flat-normalised coefficients of zeta3, pi^2 log2, log^3 2:', [sp.simplify(t) for t in tot])

# (2) squircle
mp.mp.dps = 80
b = mp.sqrt(2) - 1; r = 1/mp.sqrt(2); lam = mp.log(2); L = mp.log(1 + mp.sqrt(2)); pi2 = mp.pi**2; z3 = mp.zeta(3)
L3 = lambda v: mp.polylog(3, v)
table = {   # value -> (zeta3 part, full rhs) from the note's Lemma lem:trilog (+ D3 for b^2)
    '-x':           (sp.Rational(7, 32), 7*z3/32 - L3(r) - pi2*lam/48 + lam**3/24, L3(-r)),
    '1-x':          (1, z3 - L3(-b) - L3(r) - pi2*lam/12 + lam**3/24 + lam**2*L/8, L3(1 - r)),
    '1/(1+x)':      (sp.Rational(25, 32), 25*z3/32 - L3(b) + L3(r) + 5*pi2*lam/48 - pi2*L/6 - lam**3/48
                     - lam*L**2/4 + L**3/3, L3(2 - mp.sqrt(2))),
    '(1+x)/2':      (sp.Rational(-25, 32), -25*z3/32 - 13*L3(b) + 7*L3(-b) + 18*L3(r) + 13*pi2*lam/16
                     - 7*pi2*L/6 + 3*lam**3/16 - 9*lam**2*L/8 - 3*lam*L**2/4 + 7*L**3/3, L3((2 + mp.sqrt(2))/4)),
    '(1-x)/(1+x)':  (0, 4*L3(b) + 4*L3(-b), L3(b**2)),
    '-(1-x)/(1+x)': (sp.Rational(-11, 8), -11*z3/8 + 6*L3(b) + 6*L3(-b) + pi2*L/6 - L**3/3, L3(-b**2)),
    '2x/(1+x)':     (sp.Rational(-3, 8), -3*z3/8 - 6*L3(b) + 2*L3(-b) + 8*L3(r) + 7*pi2*lam/12 - 2*pi2*L/3
                     - lam**3/6 - lam*L**2 + 4*L**3/3, L3(2*b)),
    'x':            (0, L3(r), L3(r)),
}
Zsq = Z.subs(c, 1/sp.sqrt(2))
for arg, (zc, rhs, lhs) in table.items():
    assert abs(rhs - lhs) < mp.mpf(10)**-60, arg
    Zsq += C[arg].subs(c, 1/sp.sqrt(2))*zc
print('(2) note-table identities re-checked to 60 digits: OK')
print('    zeta(3) coefficient after reduction to {Li3(b),Li3(-b),Li3(1/sqrt2)}:', sp.radsimp(sp.nsimplify(sp.expand(Zsq))),
      '  (note: 189 sqrt2/64 - 19/16 =', sp.Rational(189, 64)*sp.sqrt(2) - sp.Rational(19, 16), ')')

# (3) 2pi/3
mp.mp.dps = 220
l2, l3 = mp.log(2), mp.log(3)
basisA = [('Li3(1/3)', L3(mp.mpf(1)/3)), ('Li3(-1/2)', L3(-mp.mpf(1)/2)), ('zeta3', mp.zeta(3)),
          ('pi2*log2', mp.pi**2*l2), ('pi2*log3', mp.pi**2*l3), ('log2^3', l2**3), ('log2^2*log3', l2**2*l3),
          ('log2*log3^2', l2*l3**2), ('log3^3', l3**3)]
vals = {'x': mp.mpf(1)/2, '-x': -mp.mpf(1)/2, '1-x': mp.mpf(1)/2, '1/(1+x)': mp.mpf(2)/3, '(1+x)/2': mp.mpf(3)/4,
        '(1-x)/(1+x)': mp.mpf(1)/3, '-(1-x)/(1+x)': -mp.mpf(1)/3, '2x/(1+x)': mp.mpf(2)/3}
Z23 = Z.subs(c, sp.Rational(1, 2))
for arg, v in vals.items():
    vec = [L3(v)] + [bv for _, bv in basisA]
    mp.mp.dps = 200
    rel = mp.pslq(vec, maxcoeff=10**8, maxsteps=10**6)
    mp.mp.dps = 220
    zc = sp.Rational(-rel[3], rel[0])
    print(f'    Li3({mp.nstr(v, 6)}) in basis A: relation {rel}, zeta3 part {zc}')
    Z23 += C[arg].subs(c, sp.Rational(1, 2))*zc
print('(3) 2pi/3: zeta(3) coefficient after reduction to basis A:', sp.nsimplify(Z23), '  (direct PSLQ: 3/16)')
