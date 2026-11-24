"""Statement check of thm:hexagon5: direct P^{(5)}(1) (own polar jets) vs the 307-monomial closed form
(own HPLs, own edge data), edge by edge.  Usage: run_hex5.py DPS [case ...]"""
import sys, time
import sympy as sp, mpmath as mp
import direct5 as D
import formula5 as FM

s3 = sp.sqrt(3)
phi = (1 + sp.sqrt(5)) / 2


def lin(M, P):
    return (M[0][0] * P[0] + M[0][1] * P[1], M[1][0] * P[0] + M[1][1] * P[1])


CASES = {
    'regular': ([(0, 1 / s3), (-sp.Rational(1, 2), 1 / (2 * s3)), (-sp.Rational(1, 2), -1 / (2 * s3))], [1, 1, 1]),
    'paper_2_half_1': ([(1, 0), (2, 2), (0, 1)], [1, 3, 1]),
    'generic_A': ([lin([[1, -2], [1, 1]], P) for P in [(1, 0), (sp.Rational(13, 10), 2), (0, 1)]], [5, 1, 2]),
    'generic_B_r4': ([(1, 0), (4, sp.Rational(16, 5)), (0, 1)], [1, 2, 3]),
    'golden': ([(1, 0), (phi, 1), (0, 1)], [2, 3, 5]),
    'generic_C': ([lin([[2, 1], [sp.Rational(-1, 3), 1]], P) for P in [(1, 0), (sp.Rational(3, 2), sp.Rational(5, 7)), (0, 1)]], [sp.Rational(3, 2), 4, 1]),
}


def run(name, dps):
    mp.mp.dps = dps
    V, w = CASES[name]
    eds = FM.edges_data(V, w)
    H = D.setup([(FM.N(a), FM.N(b)) for a, b in V], [FM.N(x) for x in w])
    print('=== %s  corners %s  weights %s' % (name, V, w))
    print('    ratios r =', [sp.nsimplify(e['r']) for e in eds])
    t0 = time.time()
    tot, per = D.all_derivs(H)
    print('    direct: P(1)=%s  P5=%s   (%.1f s, levels %s)' % (mp.nstr(tot[0], 20), mp.nstr(tot[5], dps), time.time() - t0, [p[2] for p in per]))
    t0 = time.time()
    Fsum = mp.mpf(0)
    worst = mp.mpf(0)
    for q, ed in enumerate(eds):
        if sp.simplify(ed['r'] - 1) == 0:
            Ef = FM.E_formula_r1(ed)
        else:
            Ef = FM.E_formula(ed)
        Ed = per[q][0][5]
        diff = Ef - Ed
        worst = max(worst, abs(diff))
        print('    edge %d r=%s: direct %s  formula %s  diff %s  (quad err est %s)' % (
            q, sp.nsimplify(ed['r']), mp.nstr(Ed, 30), mp.nstr(Ef, 30), mp.nstr(diff, 3), mp.nstr(per[q][1], 2)))
        Fsum += 2 * Ef
    print('    P5 formula %s   P5 direct - formula = %s   (%.1f s)' % (mp.nstr(Fsum, dps), mp.nstr(tot[5] - Fsum, 3), time.time() - t0))
    rel = abs(tot[5] - Fsum) / abs(tot[5])
    print('    agreement: %.1f digits (relative), worst edge diff %s' % (-mp.log10(rel) if rel else dps, mp.nstr(worst, 3)))
    return tot, Fsum


if __name__ == '__main__':
    dps = int(sys.argv[1])
    names = sys.argv[2:] or list(CASES)
    for nm in names:
        run(nm, dps)
        sys.stdout.flush()
