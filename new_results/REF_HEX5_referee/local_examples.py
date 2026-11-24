"""prop:order5 on fully local polygons (own code):
 1. exact P^{(5)}(1) per edge (exact_local_edge.py) vs direct polar jets (direct5.py);
 2. its zeta(5), zeta(2)zeta(3) coefficients vs the printed Z_e^{(5)}, Z_e^{(2,3)};
 3. integer relations (PSLQ) on the direct per-edge numbers alone, at 2 precisions;
 4. the printed sums C_m^{(5)}, C_m^{(2,3)} for regular 2m-gons, m = 2..10.
"""
import sys, time
import sympy as sp, mpmath as mp
import direct5 as D
import formula5 as FM
import exact_local_edge as XL
from local_law import printed, A0 as sA0, A1 as sA1, L as sL, L0 as sL0, h as sh

P5f, P23f = printed()


def Zprinted(ed):
    L = sp.sqrt(ed['L2'])
    L0 = -ed['q0'] / L
    h2 = ed['V02'] - ed['q0'] ** 2 / ed['L2']
    sub = {sA0: ed['A0'], sA1: ed['A1'], sL: L, sL0: L0, sh: sp.sqrt(h2)}
    return sp.simplify(P5f.subs(sub)), sp.simplify(P23f.subs(sub))


def coeffs(expr):
    e = sp.expand(expr)
    c5 = sp.simplify(e.coeff(sp.zeta(5)))
    c23 = sp.simplify(6 * e.coeff(sp.zeta(3)).coeff(sp.pi, 2).subs(sp.log(2), 0))
    return c5, c23


CASES = {
    'parallelogram_(3,1)_(1,2)': ([(3, 1), (1, 2)], [1, 1]),
    'parallelogram_(2,0)_(1,2)': ([(2, 0), (1, 2)], [1, 1]),
    'affine_regular_hexagon_w124': ([(1, 0), (1, 1), (0, 1)], [1, 2, 4]),
}


def run_case(name, dps):
    V, w = CASES[name]
    eds = FM.edges_data(V, w)
    m = len(V)
    mp.mp.dps = dps
    H = D.setup([(FM.N(a), FM.N(b)) for a, b in V], [FM.N(x) for x in w])
    tot, per = D.all_derivs(H)
    print('=== %s, direct P5 = %s' % (name, mp.nstr(tot[5], 40)))
    S5 = 0; S23 = 0; X5 = 0; X23 = 0; Pex = 0
    for q, ed in enumerate(eds):
        B = ed.get('B', 0)
        if m == 3:
            assert sp.simplify(ed['r'] - 1) == 0, 'not fully local'
        else:
            pass  # parallelogram: no third term
        E = XL.edge_exact(ed['A0'], ed['A1'], B, ed['L2'], ed['q0'], ed['V02'])[5]
        sub = {XL.m0: sp.log(ed['A0'] / ed['wA0']), XL.m1: sp.log(ed['A1'] / ed['wA1'])}
        sub[XL.mB] = sp.log(B / ed['wB']) if m == 3 else 0
        E = sp.expand(sp.expand_log(E.subs(sub), force=True))
        Ed = per[q][0][5]
        En = mp.mpf(sp.N(E, dps + 10))
        c5, c23 = coeffs(E)
        z5, z23 = Zprinted(ed)
        print('  edge %d: direct %s  exact %s  diff %s' % (q, mp.nstr(Ed, 30), mp.nstr(En, 30), mp.nstr(Ed - En, 3)))
        print('          [z5] exact %s   printed Z5 %s   diff %s' % (c5, z5, sp.simplify(c5 - z5)))
        print('          [z2z3] exact %s printed Z23 %s  diff %s' % (c23, z23, sp.simplify(c23 - z23)))
        S5 += 2 * z5; S23 += 2 * z23; X5 += 2 * c5; X23 += 2 * c23; Pex += 2 * E
    print('  sum: [z5] %s = %s ; [z2z3] %s = %s' % (sp.nsimplify(sp.simplify(X5)), sp.simplify(X5 - S5) == 0,
                                                   sp.nsimplify(sp.simplify(X23)), sp.simplify(X23 - S23) == 0))
    print('  P5 direct - exact = %s' % mp.nstr(tot[5] - mp.mpf(sp.N(Pex, dps + 10)), 3))
    return eds, per


def pslq_edge(name, q, dps_list=(110, 150)):
    """PSLQ of X = E_e / L on weight<=5 constants (log 2 included), from the direct numbers only."""
    V, w = CASES[name]
    eds = FM.edges_data(V, w)
    ed = eds[q]
    out = []
    for dps in dps_list:
        mp.mp.dps = dps
        H = D.setup([(FM.N(a), FM.N(b)) for a, b in V], [FM.N(x) for x in w])
        res, err, lev = D.edge_derivs(H, q)
        Xv = res[5] / mp.sqrt(FM.N(ed['L2']))
        l2 = mp.log(2); pi2 = mp.pi ** 2; z3 = mp.zeta(3); z5 = mp.zeta(5)
        if len(V) == 2:
            names = ['1', 'pi2', 'z3', 'pi4', 'z5', 'z2z3']
            basis = [1, pi2, z3, pi2 ** 2, z5, pi2 / 6 * z3]
        else:
            names = ['1', 'l2', 'pi2', 'l2^2', 'z3', 'pi2 l2', 'l2^3', 'pi4', 'z3 l2', 'pi2 l2^2', 'l2^4',
                     'z5', 'z2z3', 'pi4 l2', 'z3 l2^2', 'pi2 l2^3', 'l2^5']
            basis = [1, l2, pi2, l2 ** 2, z3, pi2 * l2, l2 ** 3, pi2 ** 2, z3 * l2, pi2 * l2 ** 2, l2 ** 4,
                     z5, pi2 / 6 * z3, pi2 ** 2 * l2, z3 * l2 ** 2, pi2 * l2 ** 3, l2 ** 5]
        rel = mp.pslq([Xv] + basis, maxcoeff=10 ** 12, maxsteps=10 ** 6, tol=mp.mpf(10) ** (-(dps - 15)))
        if rel is None:
            print('  PSLQ %s edge %d dps %d: none' % (name, q, dps)); out.append(None); continue
        resid = mp.fsum(c * v for c, v in zip(rel, [Xv] + basis))
        c = {nm: sp.Rational(-rel[i + 1], rel[0]) for i, nm in enumerate(names)}
        print('  PSLQ %s edge %d dps %d: maxcoeff %d resid %s ; [z5]=%s*L  [z2z3]=%s*L' % (
            name, q, dps, max(abs(x) for x in rel), mp.nstr(resid, 3), c['z5'], c['z2z3']))
        out.append((rel, c))
    z5, z23 = Zprinted(ed)
    L = sp.sqrt(ed['L2'])
    if out[-1] is not None:
        print('    printed Z5/L = %s, Z23/L = %s ; match: %s %s ; same relation at both precisions: %s' % (
            sp.simplify(z5 / L), sp.simplify(z23 / L), sp.simplify(out[-1][1]['z5'] - z5 / L) == 0,
            sp.simplify(out[-1][1]['z2z3'] - z23 / L) == 0, out[0] is not None and list(out[0][0]) == list(out[-1][0])))


def regular_sums():
    """numerical (50 digits): edge data of the paper's regular 2m-gon K_1 computed in floating point."""
    print('=== regular 2m-gons: sum_e Z_e (printed law) vs printed C_m')
    f5 = sp.lambdify((sA0, sA1, sL, sL0, sh), P5f, 'mpmath')
    f23 = sp.lambdify((sA0, sA1, sL, sL0, sh), P23f, 'mpmath')
    mp.mp.dps = 50
    for m in range(2, 11):
        n = [(mp.cos(mp.pi * i / m), mp.sin(mp.pi * i / m)) for i in range(1, m + 1)]
        ang = [mp.pi / 2 + mp.pi * q / m for q in range(m)]
        R = 1 / mp.fsum(abs(nn[0] * mp.cos(ang[0]) + nn[1] * mp.sin(ang[0])) for nn in n)
        V = [(R * mp.cos(a), R * mp.sin(a)) for a in ang] + [(-R * mp.cos(ang[0]), -R * mp.sin(ang[0]))]
        val = lambda x: mp.fsum(abs(nn[0] * x[0] + nn[1] * x[1]) for nn in n)
        S5 = mp.mpf(0); S23 = mp.mpf(0)
        for q in range(m):
            P0, P1 = V[q], V[q + 1]
            assert abs(val(P0) - 1) < 1e-40
            D_ = (P1[0] - P0[0], P1[1] - P0[1])
            L = mp.sqrt(D_[0] ** 2 + D_[1] ** 2)
            q0 = P0[0] * D_[0] + P0[1] * D_[1]
            L0 = -q0 / L
            h = abs(P0[0] * P1[1] - P0[1] * P1[0]) / L
            # xi0 = term of the kink line through P0, at P1; that line has normal n_j with <n_j,P0> = 0
            j0 = min(range(m), key=lambda j: abs(n[j][0] * P0[0] + n[j][1] * P0[1]))
            j1 = min(range(m), key=lambda j: abs(n[j][0] * P1[0] + n[j][1] * P1[1]))
            xi0 = abs(n[j0][0] * P1[0] + n[j0][1] * P1[1]); xi1 = abs(n[j1][0] * P0[0] + n[j1][1] * P0[1])
            S5 += 2 * f5(xi0, xi1, L, L0, h); S23 += 2 * f23(xi0, xi1, L, L0, h)
        s = mp.sin(mp.pi / (2 * m)); c = mp.cos(mp.pi / (2 * m))
        Cm5 = 2880 * m * s ** 4 * c * (65 * s ** 6 - 63 * s ** 4 + 17 * s ** 2 - 1)
        Cm23 = -960 * m * c * s ** 4 * (45 * s ** 6 - 27 * s ** 4 + 6 * s ** 2 - 1)
        print('  m=%d: sum Z5 = %s, C_m5 = %s, diff %s ; sum Z23 = %s, C_m23 = %s, diff %s' % (
            m, mp.nstr(S5, 20), mp.nstr(Cm5, 20), mp.nstr(S5 - Cm5, 3), mp.nstr(S23, 20), mp.nstr(Cm23, 20), mp.nstr(S23 - Cm23, 3)))
    print('  m=3 printed: 2835 sqrt3/32 = %s, 1395 sqrt3/32 = %s; m=2: -90 sqrt2 = %s, -210 sqrt2 = %s' % (
        mp.nstr(2835 * mp.sqrt(3) / 32, 20), mp.nstr(1395 * mp.sqrt(3) / 32, 20), mp.nstr(-90 * mp.sqrt(2), 20), mp.nstr(-210 * mp.sqrt(2), 20)))


if __name__ == '__main__':
    t0 = time.time()
    regular_sums()
    for nm in CASES:
        run_case(nm, 60)
        sys.stdout.flush()
    pslq_edge('parallelogram_(3,1)_(1,2)', 0)
    pslq_edge('parallelogram_(3,1)_(1,2)', 1)
    for q in range(3):
        pslq_edge('affine_regular_hexagon_w124', q, (160, 200))
    print('time %.0f s' % (time.time() - t0))
