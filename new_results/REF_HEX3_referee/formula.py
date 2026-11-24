"""Referee's transcription of the K1e coefficient table (formula_table.txt, text only) and of the
printed Proposition prop:nonlocal / Proposition prop:edgelaw formulas of main.tex.
Edge data are computed here from scratch from the corners and weights.
"""
import re, sympy as sp, mpmath as mp

TABLE = '<project>/research/08_theory/K1e_hexagon_closed/formula_table.txt'
B, r, q0, L2, V02, lA0, lA1, lB = sp.symbols('B r q0 L2 V02 lA0 lA1 lB')
SYMS = dict(B=B, r=r, q0=q0, L2=L2, V02=V02, lA0=lA0, lA1=lA1, lB=lB)

def parse_table():
    blocks = {}
    cur = None
    for line in open(TABLE):
        m = re.match(r'===== coefficient of (.+?)\s+\(denominator', line)
        if m:
            cur = m.group(1).strip(); blocks[cur] = sp.Integer(0); continue
        m = re.match(r'\s*\[(.+?)\] \* \( (.+) \)\s*$', line)
        if m and cur is not None:
            mono = sp.sympify(m.group(1), locals=SYMS)
            expr = sp.sympify(m.group(2), locals=SYMS)
            blocks[cur] += mono * expr
    return blocks

BLOCKS = parse_table()

def edges_data(V, w):
    """Exact (sympy) edge data for edges V0->V1, V1->V2, V2->-V0 of the hexagon with weights w."""
    V = [sp.Matrix([sp.nsimplify(a), sp.nsimplify(b)]) for a, b in V]
    cr = lambda P, Q: P[0] * Q[1] - P[1] * Q[0]
    M = sp.Matrix(3, 3, lambda q, j: abs(cr(V[j], V[q])))
    c = list(M.LUsolve(sp.Matrix([1, 1, 1])))
    w = [sp.nsimplify(x) for x in w]
    ell = lambda j, X: c[j] * abs(cr(V[j], X))
    out = []
    corners = [(V[0], 0), (V[1], 1), (V[2], 2), (-V[0], 0)]
    for q in range(3):
        (P0, j0), (P1, j1) = corners[q], corners[q + 1]
        jB = 3 - j0 - j1
        A0 = sp.nsimplify(ell(j0, P1)); A1 = sp.nsimplify(ell(j1, P0)); Bv = sp.nsimplify(ell(jB, P0))
        rv = sp.nsimplify(ell(jB, P1) / Bv)
        assert sp.simplify(A1 + Bv - 1) == 0 and sp.simplify(A0 - (1 - Bv * rv)) == 0
        D = P1 - P0
        out.append(dict(B=Bv, r=rv, q0=sp.nsimplify(P0.dot(D)), L2=sp.nsimplify(D.dot(D)),
                        V02=sp.nsimplify(P0.dot(P0)), A0=A0, A1=A1,
                        wA0=w[j0], wA1=w[j1], wB=w[jB]))
    return out

def mp_of(e):
    return mp.mpf(sp.N(e, mp.mp.dps + 20)) if not isinstance(e, mp.mpf) else e

def E_formula(ed, sign_convention='ell'):
    """E_e from the table.  'ell': the label ln(1-x) is read as ell=-log r (planner correction, and main.tex);
    'literal': ln(1-x)=log r as the table labels literally say."""
    Bv, rv = mp_of(ed['B']), mp_of(ed['r'])
    x = 1 - rv
    ell = -mp.log(rv) if sign_convention == 'ell' else mp.log(rv)
    vals = {B: Bv, r: rv, q0: mp_of(ed['q0']), L2: mp_of(ed['L2']), V02: mp_of(ed['V02']),
            lA0: mp.log(mp_of(ed['A0']) / mp_of(ed['wA0'])), lA1: mp.log(mp_of(ed['A1']) / mp_of(ed['wA1'])),
            lB: mp.log(Bv / mp_of(ed['wB']))}
    atoms = {
        'zeta3': mp.zeta(3), 'Li3(x)': mp.polylog(3, x), 'Li3(x/(x-1))': mp.polylog(3, x / (x - 1)),
        'Li2(x)*ln(1-x)': mp.polylog(2, x) * ell, 'zeta2*ln(1-x)': mp.zeta(2) * ell, 'zeta2': mp.zeta(2),
        'Li2(x)': mp.polylog(2, x), 'ln^3(1-x)': ell ** 3, 'ln^2(1-x)': ell ** 2, 'ln(1-x)': ell, '1': mp.mpf(1)}
    tot = 0
    for name, expr in BLOCKS.items():
        f = sp.lambdify(list(vals.keys()), expr, 'mpmath')
        tot += f(*vals.values()) * atoms[name]
    return tot / vals[L2] ** mp.mpf(2.5)

# ---- printed formulas of main.tex ----
def Ze_printed(ed):
    """prop:edgelaw: Z_e = 6 xi0 xi1 h^2/L (3 (xi0 L0 + xi1 L1)/L - 2)."""
    L = sp.sqrt(ed['L2']); h2 = (ed['L2'] * ed['V02'] - ed['q0'] ** 2) / ed['L2']
    L0 = -ed['q0'] / L; L1 = L - L0
    xi0, xi1 = ed['A0'], ed['A1']
    return sp.nsimplify(sp.radsimp(6 * xi0 * xi1 * h2 / L * (3 * (xi0 * L0 + xi1 * L1) / L - 2)))

def Phi_printed(ed):
    L = sp.sqrt(ed['L2']); h2 = (ed['L2'] * ed['V02'] - ed['q0'] ** 2) / ed['L2']
    Bv, rv, q = ed['B'], ed['r'], ed['q0']
    A0, A1, kap = ed['A0'], ed['A1'], ed['r'] - 1
    gam = Bv * h2 / L ** 3
    PA = -6 * gam * (3 * A0 * kap * q - L ** 2 * (2 * A0 + rv * A1))
    PL = -18 * gam * (kap * (-A1 ** 2 + A1 * Bv * (2 * kap - 1) + 2 * Bv ** 2 * kap) * q + L ** 2 * (A0 * Bv + rv * A1 ** 2))
    return sp.radsimp(PA), sp.radsimp(PL)

def Ne_printed(ed):
    PA, PL = Phi_printed(ed)
    rv = ed['r']
    if rv == 2:
        zA, zL = sp.Rational(-5, 8), sp.Rational(-3, 4)
    elif rv == sp.Rational(1, 2):
        zA, zL = sp.Integer(1), sp.Rational(7, 8)
    elif rv == 1:
        return sp.Integer(0)
    else:
        raise ValueError('prop:nonlocal only for r in {2,1/2,1}')
    return sp.radsimp(PA * zA + PL * zL)
