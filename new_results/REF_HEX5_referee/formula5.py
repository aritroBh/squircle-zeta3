"""Evaluate the closed form of thm:hexagon5 from the notes' 307-monomial table
(K1g_hexagon_order5_assembly/coef_acc_raw.pkl, read only: a dict
 (atom-monomial, x-power, (z2,z3,z5)-powers) -> coefficient in A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB),
with the referee's own HPLs (hpl.py) and the referee's own edge data.
E_e = 120 * sum coef * atoms * x^xe * zeta's,  x = -kap = 1-r;  P^{(5)}(1) = 2 sum_e E_e.
"""
import pickle, os, time
import sympy as sp, mpmath as mp
import hpl

TAB = '<project>/research/08_theory/K1g_hexagon_order5_assembly/coef_acc_raw.pkl'
A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
SYMS = [A0, A1, B, kap, L2, q0, V02, lA0, lA1, lB]
_F = None


def load():
    global _F
    if _F is None:
        t0 = time.time()
        coef = pickle.load(open(TAB, 'rb'))
        _F = []
        for key, co in coef.items():
            f = sp.lambdify(SYMS, co, modules='mpmath')
            _F.append((key, f))
        print('# table: %d monomials, lambdified in %.1f s' % (len(_F), time.time() - t0))
    return _F


def edges_data(V, w):
    """Exact edge data for edges V0->V1, V1->V2, ..., V_{m-1}->-V0 (hexagon: m = 3).
    Term of kink line j at x: c_j |cross(V_j, x)| (sum = 1 on the boundary at p = 1)."""
    V = [sp.Matrix([sp.nsimplify(a), sp.nsimplify(b)]) for a, b in V]
    m = len(V)
    cr = lambda P, Q: P[0] * Q[1] - P[1] * Q[0]
    M = sp.Matrix(m, m, lambda q, j: abs(cr(V[j], V[q])))
    c = list(M.LUsolve(sp.Matrix([1] * m)))
    c = [sp.nsimplify(sp.simplify(x)) for x in c]
    w = [sp.nsimplify(x) for x in w]
    ell = lambda j, X: sp.simplify(c[j] * abs(cr(V[j], X)))
    corners = [(V[q], q) for q in range(m)] + [(-V[0], 0)]
    out = []
    for q in range(m):
        (P0, j0), (P1, j1) = corners[q], corners[q + 1]
        others = [j for j in range(m) if j not in (j0, j1)]
        a0 = ell(j0, P1)   # A0 t : term of the line through V0, value at V1
        a1 = ell(j1, P0)   # A1 (1-t)
        D = P1 - P0
        d = dict(A0=a0, A1=a1, q0=sp.simplify(P0.dot(D)), L2=sp.simplify(D.dot(D)),
                 V02=sp.simplify(P0.dot(P0)), wA0=w[j0], wA1=w[j1], c=c, j0=j0, j1=j1)
        if m == 3:
            jB = others[0]
            Bv = ell(jB, P0)
            rv = sp.simplify(ell(jB, P1) / Bv)
            assert sp.simplify(a1 + Bv - 1) == 0 and sp.simplify(a0 - (1 - Bv * rv)) == 0
            d.update(B=Bv, r=rv, wB=w[jB])
        out.append(d)
    return out


def N(e, dps=None):
    return mp.mpf(sp.N(e, (dps or mp.mp.dps) + 20))


def E_formula(ed, x_override=None, guard=40):
    with mp.workdps(mp.mp.dps + guard):
        v = _E_formula(ed, x_override)
    return +v


def _E_formula(ed, x_override=None):
    """E_e of the closed form.  If x_override is given, evaluate at x = x_override (kap = -x) with
    A0 = 1 - B(1+kap) and all other data fixed (used for r = 1 edges, where the monomials have poles)."""
    F = load()
    Bv = N(ed['B'])
    if x_override is None:
        rv = N(ed['r'])
        xv = 1 - rv
        A0v = N(ed['A0'])
    else:
        xv = mp.mpf(x_override)
        A0v = 1 - Bv * (1 - xv)
    kv = -xv
    args = [A0v, N(ed['A1']), Bv, kv, N(ed['L2']), N(ed['q0']), N(ed['V02']),
            mp.log(N(ed['A0']) / N(ed['wA0'])), mp.log(N(ed['A1']) / N(ed['wA1'])), mp.log(Bv / N(ed['wB']))]
    at = hpl.atoms(xv)
    z = (mp.zeta(2), mp.zeta(3), mp.zeta(5))
    tot = mp.mpf(0)
    for (atm, xe, ze), f in F:
        mv = xv ** xe * z[0] ** ze[0] * z[1] ** ze[1] * z[2] ** ze[2]
        for nm, e in atm:
            mv *= at[nm] ** e
        tot += f(*args) * mv
    return 120 * tot


def E_formula_r1(ed, eps_exp=None):
    """r = 1 edge: average of the closed form at x = +eps and -eps (error O(eps^2))."""
    dps = mp.mp.dps
    if eps_exp is None:
        eps_exp = dps // 2 + 5
    eps = mp.mpf(10) ** (-eps_exp)
    with mp.workdps(dps + 7 * eps_exp + 20):
        v = (E_formula(ed, eps) + E_formula(ed, -eps)) / 2
    return +v
