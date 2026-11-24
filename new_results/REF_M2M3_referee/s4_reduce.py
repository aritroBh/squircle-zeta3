"""Referee step 2c: exact reduction of the Li_2 / Li_3 values produced by s3_weight3.py, using only
instances (inside their real domains) of classical functional equations, and comparison with
(T1), (T2), T = (pi^2+6 lambda^2)/24, and the note's Lemma dilog / Lemma trilog.

Usage: python s4_reduce.py note   -> Kummer (K) only at the four points used in the note
       python s4_reduce.py full   -> all (K) and Abel instances whose arguments lie in the argument set
       python s4_reduce.py noK    -> no Kummer at all (to see whether K is really needed)
"""
import sys
import itertools
import sympy as sp
import mpmath as mp
from sympy.polys.matrices import DomainMatrix
import s3_weight3 as W
from q2 import Q2, A, B, R2, ONE, ZERO, logS, lam, Lg

MODE = sys.argv[1] if len(sys.argv) > 1 else 'note'
mp.mp.dps = 80
pi2 = sp.pi**2
z3 = W.z3
Li = W.Li
lg = logS
s2 = sp.sqrt(2)
K2 = sp.QQ.algebraic_field(s2)

Tx = W.integrate(W.T_terms)
W3x = W.integrate(W.W3_terms)
W3px = W.integrate(W.W3p_terms)

R2inv = ONE/R2
HALF = Q2(sp.Rational(1, 2).p, 1) / 2 if False else Q2(1)/2
basis = {2: [B], 3: [B, -B, R2inv]}

# ---------------- argument set: arguments that occur, basis, specials; closed under z->1-z, 1/z
S = {u for (n, u) in W.REG} | {B, -B, R2inv, ONE, -ONE, HALF}


def anh(u):
    out = [u]
    if u != ONE and u != ZERO:
        out += [ONE - u, ONE/u, ONE/(ONE - u), u/(u - ONE), (u - ONE)/u]
    return [w for w in out if w != ZERO and w != ONE] + ([ONE] if u == ONE else [])


changed = True
while changed:
    changed = False
    for u in list(S):
        for w in anh(u):
            if w not in S:
                S.add(w); changed = True
print('argument set size', len(S))


def isS(u):
    try:
        logS(u); return True
    except AssertionError:
        return False


special = {
    2: {ONE: pi2/6, -ONE: -pi2/12, HALF: pi2/12 - lam**2/2},
    3: {ONE: z3, -ONE: -sp.Rational(3, 4)*z3, HALF: sp.Rational(7, 8)*z3 - pi2/12*lam + lam**3/6},
}

rels = {2: [], 3: []}
names = {2: [], 3: []}


def add(n, expr, name, check_points):
    rels[n].append(sp.expand(expr)); names[n].append(name)


for n in (2, 3):
    for u, val in special[n].items():
        add(n, Li(n, u) - val, f'special Li{n}({u})', None)

for z in S:
    # (R) 0<z<1
    if ZERO < z < ONE and (ONE - z) in S:
        add(2, Li(2, z) + Li(2, ONE - z) - pi2/6 + lg(z)*lg(ONE - z), f'R({z})', None)
    # (I) z<0
    if z < ZERO and (ONE/z) in S:
        add(2, Li(2, z) + Li(2, ONE/z) + pi2/6 + lg(-z)**2/2, f'I({z})', None)
        add(3, Li(3, z) - Li(3, ONE/z) + pi2/6*lg(-z) + lg(-z)**3/6, f'I3({z})', None)
    # (La) z<1
    if z < ONE and z/(z - ONE) in S:
        add(2, Li(2, z) + Li(2, z/(z - ONE)) + lg(ONE - z)**2/2, f'La({z})', None)
    # (T3) 0<z<1
    if ZERO < z < ONE and (ONE - z) in S and (ONE - ONE/z) in S:
        add(3, Li(3, z) + Li(3, ONE - z) + Li(3, ONE - ONE/z)
            - (z3 + lg(z)**3/6 + pi2/6*lg(z) - lg(z)**2*lg(ONE - z)/2), f'T3({z})', None)
    # (D),(D3) |w|<=1
    if abs(z) <= ONE and z*z in S and -z in S:
        add(2, Li(2, z*z) - 2*Li(2, z) - 2*Li(2, -z), f'D({z})', None)
        add(3, Li(3, z*z) - 4*Li(3, z) - 4*Li(3, -z), f'D3({z})', None)

# Abel five-term (A): x<1, y<1, x+y<1
for x, y in itertools.product(S, S):
    if x < ONE and y < ONE and x + y < ONE and x != ZERO and y != ZERO:
        args = [x/(ONE - y), y/(ONE - x), x*y/((ONE - x)*(ONE - y))]
        if all(u in S or u == ZERO for u in args):
            add(2, Li(2, args[0]) + Li(2, args[1]) - Li(2, args[2]) - Li(2, x) - Li(2, y)
                - lg(ONE - x)*lg(ONE - y), f'A({x},{y})', None)


def kummer(x, y):
    args = [x*(ONE - y)*(ONE - y)/(y*(ONE - x)*(ONE - x)), x*y, x/y, x*(ONE - y)/(y*(ONE - x)),
            x*(ONE - y)/(x - ONE), (ONE - y)/(ONE - x), (ONE - y)/(y*(x - ONE)), x, y]
    cf = [1, 1, 1, -2, -2, -2, -2, -2, -2]
    lhs = sum(c*Li(3, u) for c, u in zip(cf, args)) + 2*z3
    rhs = lg(y)**2*lg((ONE - y)/(ONE - x)) - pi2/3*lg(y) - lg(y)**3/3
    return args, lhs - rhs


NOTE_K = [(B, R2inv), (Q2(2) - R2, R2inv), (R2inv, 2*B), (-ONE, B*B)]
if MODE == 'note':
    Kpts = NOTE_K
elif MODE == 'full':
    Kpts = [(x, y) for x, y in itertools.product(S, S) if ZERO < y < ONE and x < y and x != ZERO]
else:
    Kpts = []
nK = 0
for x, y in Kpts:
    assert ZERO < y < ONE and x < y        # domain D of the real form
    args, rel = kummer(x, y)
    if all(u in S or u == ZERO for u in args):
        add(3, rel, f'K({x},{y})', None); nK += 1
    else:
        assert MODE != 'note', ('note K point leaves argument set', x, y)
print('relations: weight 2:', len(rels[2]), ' weight 3:', len(rels[3]), ' (K instances:', nK, ')')

# ---------------- numerical sanity check of every relation instance (catches misquoted equations)
def numval(expr):
    return W.numval(expr)


worst = {2: 0, 3: 0}
for n in (2, 3):
    for r, nm in zip(rels[n], names[n]):
        d = abs(numval(r))
        worst[n] = max(worst[n], d)
        if d > mp.mpf(10)**-50:
            print('  FAILED numerically:', nm, mp.nstr(d, 3))
print('max |relation| numerically (dps 60): w2', mp.nstr(worst[2], 3), ' w3', mp.nstr(worst[3], 3))


# ---------------- exact elimination
def solve(n):
    lisyms = {u: Li(n, u) for u in S if u != ZERO}
    bas = [lisyms[u] for u in basis[n]]
    nonbas = [s for u, s in lisyms.items() if s not in bas]
    prods = []
    rows = []
    gens = nonbas + bas + [lam, Lg, sp.pi, z3]
    polys = [sp.Poly(r, *gens, domain=K2) for r in rels[n]]
    colkey = {}
    cols = []
    for s_ in nonbas + bas:
        mon = tuple(1 if g == s_ else 0 for g in gens)
        colkey[mon] = len(cols); cols.append(s_)
    for P in polys:
        for mon in P.monoms():
            if mon not in colkey:
                colkey[mon] = len(cols)
                cols.append(sp.Mul(*[g**e for g, e in zip(gens, mon)]))
    M = [[K2.zero]*len(cols) for _ in polys]
    for i, P in enumerate(polys):
        for mon, cf in P.terms():
            M[i][colkey[mon]] = cf
    DM = DomainMatrix(M, (len(polys), len(cols)), K2)
    R, piv = DM.rref()
    R = sp.Matrix([[K2.to_sympy(e) if hasattr(e, 'to_list') else sp.sympify(e) for e in row]
                   for row in R.to_list()])
    nnb = len(nonbas)
    sol = {}
    for i, pc in enumerate(piv):
        if pc < nnb:
            row = R.row(i)
            if any(row[j] != 0 for j in range(nnb) if j != pc):
                continue
            sol[cols[pc]] = sp.expand(-sum(row[j]*cols[j] for j in range(nnb, len(cols))))
        elif pc < nnb + len(bas):
            print('  WARNING: relation among basis elements found, pivot', cols[pc])
    return sol, cols


sol2, _ = solve(2)
sol3, _ = solve(3)


def reduce(expr):
    expr = sp.expand(expr)
    need = [s_ for s_ in expr.free_symbols if str(s_).startswith('Li')]
    bas = {Li(2, B), Li(3, B), Li(3, -B), Li(3, R2inv)}
    missing = [s_ for s_ in need if s_ not in bas and s_ not in sol2 and s_ not in sol3]
    out = sp.expand(expr.xreplace({**sol2, **sol3}))
    return out, missing


Li2b, Li3b, Li3mb, Li3r = Li(2, B), Li(3, B), Li(3, -B), Li(3, R2inv)
L = Lg
T_target = (pi2 + 6*lam**2)/24
T1 = (T_target + (189*s2/64 - sp.Rational(19, 16))*z3 + (11 + 15*s2)*Li3b + (11 - 15*s2)*Li3mb
      - 27*s2*Li3r + (12*L - 7*s2/2*lam)*Li2b + (15*s2/8 - sp.Rational(7, 12))*pi2*L
      - (sp.Rational(1, 24) + 65*s2/48)*pi2*lam + (sp.Rational(1, 6) + 9*s2/16)*lam**3
      + 11*s2/8*lam**2*L - (1 + 7*s2/4)*lam*L**2 + (sp.Rational(25, 6) - 7*s2/3)*L**3)
T2 = ((7*s2/32 - sp.Rational(3, 16))*z3 + (1 + s2)*Li3b + (1 - s2)*Li3mb - 2*s2*Li3r
      + (2*L - s2*lam)*Li2b + (s2/12 - sp.Rational(1, 8))*pi2*L - (1 + s2)/24*pi2*lam
      + (1 + s2)/24*lam**3 + s2/4*lam**2*L - (1 + s2)/2*lam*L**2 + (5 - s2)/6*L**3)

for name, raw, target in [('T', Tx, T_target), ('W3 = M2 - T  vs (T1) - T', W3x, T1 - T_target),
                          ("W3' vs (T2)", W3px, T2)]:
    red, missing = reduce(raw)
    diff = sp.simplify(sp.expand(sp.radsimp(red - target)))
    print(f'{name}: unreduced symbols {missing}; exact difference = {diff}')

# ---------------- the note's lemma lines, compared with my own reductions
lemma3 = {
    -A: Li3mb - pi2/6*L - L**3/6,
    ONE - R2inv: z3 - Li3mb - Li3r - pi2/12*lam + lam**3/24 + lam**2*L/8,
    -R2inv: sp.Rational(7, 32)*z3 - Li3r - pi2/48*lam + lam**3/24,
    Q2(2) - R2: sp.Rational(25, 32)*z3 - Li3b + Li3r + 5*pi2/48*lam - pi2/6*L - lam**3/48 - lam*L**2/4 + L**3/3,
    -B*B: -sp.Rational(11, 8)*z3 + 6*Li3b + 6*Li3mb + pi2/6*L - L**3/3,
    2*B: -sp.Rational(3, 8)*z3 - 6*Li3b + 2*Li3mb + 8*Li3r + 7*pi2/12*lam - 2*pi2/3*L - lam**3/6 - lam*L**2 + 4*L**3/3,
    -B/2: sp.Rational(11, 8)*z3 + 2*Li3b - 6*Li3mb - 8*Li3r - 5*pi2/12*lam + pi2/2*L + lam**3/3 + lam**2*L/2 - lam*L**2/2 - L**3/2,
    -2*A: sp.Rational(11, 8)*z3 + 2*Li3b - 6*Li3mb - 8*Li3r - 7*pi2/12*lam + pi2/3*L + lam**3/6 - lam*L**2 - 2*L**3/3,
    (Q2(2) - R2)/4: sp.Rational(101, 32)*z3 + 7*Li3b - 13*Li3mb - 18*Li3r - 17*pi2/16*lam + 7*pi2/6*L + sp.Rational(15, 16)*lam**3 + sp.Rational(9, 8)*lam**2*L - sp.Rational(3, 4)*lam*L**2 - sp.Rational(4, 3)*L**3,
    (Q2(2) + R2)/4: -sp.Rational(25, 32)*z3 - 13*Li3b + 7*Li3mb + 18*Li3r + 13*pi2/16*lam - 7*pi2/6*L + sp.Rational(3, 16)*lam**3 - sp.Rational(9, 8)*lam**2*L - sp.Rational(3, 4)*lam*L**2 + sp.Rational(7, 3)*L**3,
}
lemma2 = {
    -B: Li2b - pi2/8 + L**2/2,
    B*B: 4*Li2b - pi2/4 + L**2,
    Q2(2) - R2: -Li2b + pi2/6 + lam*L/2 - L**2,
    -R2: Li2b - pi2/6 - lam*L/2 + L**2/2,
    -B*B: 6*Li2b - sp.Rational(11, 24)*pi2 + 2*L**2,
    (Q2(2) + R2)/4: 6*Li2b - sp.Rational(7, 24)*pi2 - sp.Rational(9, 8)*lam**2 - sp.Rational(3, 2)*lam*L + sp.Rational(7, 2)*L**2,
}
for n, lem, sol in [(2, lemma2, sol2), (3, lemma3, sol3)]:
    for u, rhs in lem.items():
        sym = Li(n, u)
        num = abs(mp.polylog(n, u.mp()) - numval(rhs))
        if sym in sol:
            d = sp.simplify(sp.expand(sp.radsimp(sol[sym] - rhs)))
            print(f'  note Lemma Li{n}({u}): own exact reduction - lemma = {d};  numeric |diff| = {mp.nstr(num, 3)}')
        else:
            print(f'  note Lemma Li{n}({u}): not reduced by my system; numeric |diff| = {mp.nstr(num, 3)}')
