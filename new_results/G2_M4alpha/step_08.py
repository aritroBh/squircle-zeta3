"""step_08.py -- B4 zeta parts at weights 1..3: z3 part = -Psi(c) (B3's z2 part), and the z2-part fit.

z3 part: [z3]mu_{4,n} = -[z2]mu_{3,n} (step_05, all n<=520). B3's z2 part Psi(c) is the HPL
combination in G1's hpl_fit_k3_(1,0,0).json. Checked: the series of -Psi agrees with the exact
z3 target through c^1040 (EXACT-1040). (PROVED for all n in step_11 by the symbolic formulas.)
z2 part: fitted in the even-sector basis {c^e (H_w +/- H_{-w})}, w of weight <= 3, e in {-2..1},
by exact rational linear algebra on the coefficients 0..1040 (521 even equations, 79 unknowns);
accepted only if every equation is satisfied. Writes hpl_fit_k4_(1,0,0).json.
"""
import json
from fractions import Fraction as Fr
import g2_fit as gf

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
ORDER = 1040


def load_g1_terms(name):
    out = []
    for nm, v in json.load(open(G1 + name)):
        if nm == '1':
            out.append(('1', (), 0, Fr(v))); continue
        if nm == 'c':
            out.append(('c', (), 1, Fr(v))); continue
        e = int(nm.split(' ', 1)[0][2:])
        w = tuple(int(t) for t in nm.split(' ', 1)[1][2:-1].split(','))
        out.append(('hpl', w, e, Fr(v)))
    return out


def combo_series(terms, order):
    T = [Fr(0)]*(order + 1)
    for kind, w, e, q in terms:
        if kind == '1':
            T[0] += q
        elif kind == 'c':
            T[1] += q
        else:
            s = gf.hpl_series(w, order)
            for j in range(max(0, e), order + 1):
                m = j - e
                if s[m]:
                    T[j] += q*s[m]
    return T


if __name__ == '__main__':
    # ---- z3 part: -Psi(c)
    Psi = load_g1_terms('hpl_fit_k3_(1,0,0).json')
    sPsi = combo_series(Psi, ORDER)
    T3 = gf.target_series('(0, 1, 0)', ORDER)
    ok = all(sPsi[m] == -T3[m] for m in range(ORDER + 1))
    print('z3 part: -Psi(c) == target series through c^1040 (all coefficients, even and odd):', ok)
    if not ok:
        bad = [m for m in range(ORDER + 1) if sPsi[m] != -T3[m]][:6]
        print('  first mismatches at m =', bad)

    # ---- z2 part fit
    T2 = gf.target_series('(1, 0, 0)', ORDER)
    basis = gf.build_even_basis(3, [-3, -2, -1, 0, 1, 2], ORDER)
    terms, rankA, rank = gf.solve_exact(basis, T2, ORDER)
    if terms is None:
        print('NO SOLUTION in weight<=3, e in {-2..1} basis')
    else:
        print(f'z2 part: SOLUTION with {len(terms)} terms (rank {rankA} of {len(basis)}):')
        for nm, v in terms:
            print(f'   {v}  *  {nm}')
        json.dump([(nm, str(v)) for nm, v in terms], open('hpl_fit_k4_(1,0,0).json', 'w'))
        print('saved hpl_fit_k4_(1,0,0).json')
