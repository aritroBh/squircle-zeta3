"""step_05.py -- EXACT symbolic formula for mu_{4,n} (k=4), extending mu3_symbolic.py to weight-4 zeta content.

  mu_{4,n} = -(1/(n(2n-1))) sum_{(j,i,l)} e_{j,i,l}(n) (n)_j/(2n)_j D_{i,l}(n+j, n),
an exact polynomial in h1..h4 = H^(r)_{n-1}, g1..g4 = H^(r)_{2n-1}, z2, z3, z4 with Q(n) coefficients.
New at k=4: the third polygamma brings 6 z4 - 6 H^(4); monomials z2^2 and z4 occur (both ~ pi^4).
No z5 and no z2 z3 (Bell polynomials of degree <= 4 in the cumulants kappa_m, m <= 4).

Checks:
 (1) formula evaluation == flat_exact.py's independent direct computation, all n <= 400, exact;
 (2) [z4] = (n+3)/((2n+3)(4n^2-1)); [z2^2] = -n/((2n+3)(4n^2-1)); [z3]mu_4 = -[z2]mu_3 (all n).
Writes mu4_formula.pkl (symbolic formula, for the all-n proof) -- the table itself is
flat_exact_N400_K4.json written by the reference flat_exact.py run (see flat_exact_N400_K4.out).
"""
import sys, json, pickle, time
from fractions import Fraction as Fr
from math import comb
import sympy as sp
from sympy import QQ
from sympy.polys.rings import ring

sys.path.insert(0, "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4")
from zsum import K, nn

KMAX = 4
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
u, y, a, b, n = sp.symbols('u y a b n')

t0 = time.time()
A = -(n - 1)*a + n*b
l = b - a
logE = u*A + (1 - 2*n)*sp.log(1 + y*(sp.exp(u*l) - 1))
E_ser = sp.series(sp.exp(sp.series(logE, u, 0, KMAX + 1).removeO()), u, 0, KMAX + 1).removeO()
EK = {k: sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b) for k in range(KMAX + 1)}
print(f"symbolic E-series done ({time.time()-t0:.0f}s); terms in [u^4]E: {len(EK[4].terms())}", flush=True)

R, h1, h2, h3, h4, g1, g2, g3, g4, z2, z3, z4 = ring('h1,h2,h3,h4,g1,g2,g3,g4,z2,z3,z4', K)


def toK(expr):
    P = sp.Poly(sp.expand(expr), n)
    out = K(0)
    for (e,), c in P.terms():
        out += K(QQ(int(sp.Rational(c).p), int(sp.Rational(c).q)))*nn**e
    return out


def psiR(m, which, j):
    H = {1: h1, 2: h2, 3: h3, 4: h4} if which == 'n' else {1: g1, 2: g2, 3: g3, 4: g4}
    base = nn if which == 'n' else 2*nn
    r = m + 1
    harm = H[r] + sum((K(1)/(base + i)**r for i in range(j)), K(0))
    if m == 0:
        return harm
    if m == 1:
        return z2 - harm
    if m == 2:
        return -2*z3 + 2*harm
    if m == 3:
        return 6*z4 - 6*harm
    raise ValueError


def moments(j):
    kap = {}
    for r in range(0, KMAX + 1):
        for m in range(0, KMAX + 1 - r):
            if r + m == 0:
                continue
            if m == 0:
                kap[(r, m)] = psiR(r - 1, 'n', j) - psiR(r - 1, '2n', j)
            elif r == 0:
                kap[(r, m)] = psiR(m - 1, 'n', 0) - psiR(m - 1, '2n', j)
            else:
                kap[(r, m)] = -psiR(r + m - 1, '2n', j)
    D = {(0, 0): R(1)}
    for tot in range(1, KMAX + 1):
        for I in range(tot, -1, -1):
            L = tot - I
            if I > 0:
                v = R(0)
                for r in range(1, I + 1):
                    for m in range(0, L + 1):
                        v += comb(I - 1, r - 1)*comb(L, m)*kap[(r, m)]*D[(I - r, L - m)]
            else:
                v = R(0)
                for m in range(1, L + 1):
                    v += comb(L - 1, m - 1)*kap[(0, m)]*D[(0, L - m)]
            D[(I, L)] = v
    return D


def mu_formula(k):
    tot = R(0)
    Dc = {}
    for (j, i, ll), co in EK[k].terms():
        if j not in Dc:
            Dc[j] = moments(j)
        pr = K(1)
        for t in range(j):
            pr *= (nn + t)/(2*nn + t)
        tot += toK(co)*pr*Dc[j][(i, ll)]
    return tot*(-K(1)/(nn*(2*nn - 1)))


def Hf(r, m):
    return sum((Fr(1, i**r) for i in range(1, m + 1)), Fr(0))


def evalK(c, nv):
    return Fr(int(c.numer(nv)), 1)/Fr(int(c.denom(nv)), 1)


def zkey(mon):
    e2, e3, e4 = mon[8], mon[9], mon[10]
    assert e3 == 0 or (e2 == 0 and e4 == 0)
    assert not (e2 and e4)
    assert e2 <= 2 and e3 <= 1 and e4 <= 1
    if e2 == 1: return 'z2'
    if e2 == 2: return 'z2^2'
    if e3: return 'z3'
    if e4: return 'z4'
    return '1'


def eval_formula(F, nv, hv):
    out = {'1': Fr(0), 'z2': Fr(0), 'z3': Fr(0), 'z4': Fr(0), 'z2^2': Fr(0)}
    for mon, cf in F.terms():
        v = evalK(cf, nv)
        for x, e in zip(hv, mon[:8]):
            if e:
                v *= x**e
        out[zkey(mon)] += v
    return out


if __name__ == '__main__':
    F4 = mu_formula(4)
    print(f"mu_{{4,n}}: {len(F4.terms())} monomials ({time.time()-t0:.0f}s)", flush=True)
    pickle.dump([(mon, str(cf.as_expr())) for mon, cf in F4.terms()], open('mu4_formula.pkl', 'wb'))

    keymap = {'1': '(0, 0, 0)', 'z2': '(1, 0, 0)', 'z3': '(0, 1, 0)', 'z4': '(0, 0, 1)', 'z2^2': '(2, 0, 0)'}
    # (1) full cross-check vs flat_exact.py direct computation, all available n
    import os
    reffile = 'flat_exact_N520_K4.json' if os.path.exists('flat_exact_N520_K4.json') else 'flat_exact_N400_K4.json'
    print('cross-checking against', reffile)
    ref = json.load(open(reffile))['4']
    NREF = max(int(x) for x in ref)
    hn = {r: Fr(0) for r in range(1, 5)}
    hg = {r: Fr(0) for r in range(1, 5)}
    ok = True
    tab = {}
    for nv in range(1, NMAX + 1):
        for r in range(1, 5):
            hn[r] += Fr(1, nv**r)
        for m in range(2*nv - 1, 2*nv + 1):
            for r in range(1, 5):
                hg[r] += Fr(1, m**r)
        hv = [hn[r] - Fr(1, nv**r) for r in range(1, 5)] + [hg[r] - Fr(1, (2*nv)**r) for r in range(1, 5)]
        fv = eval_formula(F4, nv, hv)
        if nv <= NREF:
            gv = ref[str(nv)]
            for kk, gkey in keymap.items():
                if fv[kk] != Fr(gv.get(gkey, '0')):
                    ok = False; print('MISMATCH formula vs flat_exact', nv, kk)
        tab[nv] = fv
        if nv % 100 == 0:
            print(f'  n={nv} checked ({time.time()-t0:.0f}s)', flush=True)
    print('(1) formula == flat_exact.py direct, all n<=', NREF, ':', ok, flush=True)
    json.dump({'4': {str(nv): {keymap[kk]: str(v) for kk, v in tab[nv].items() if v} for nv in tab}},
              open(f'flat_exact_N{NMAX}_K4_formula.json', 'w'))
    print('saved', f'flat_exact_N{NMAX}_K4_formula.json')

    # (2) patterns
    G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
    mu3 = json.load(open(G1 + 'flat_exact_N320_K3.json'))
    p_ok = q_ok = s_ok = True
    for nv in range(1, NMAX + 1):
        D = (2*nv + 3)*(4*nv*nv - 1)
        if tab[nv]['z4'] != Fr(nv + 3, D): p_ok = False
        if tab[nv]['z2^2'] != Fr(-nv, D): q_ok = False
        if nv <= 320:
            if tab[nv]['z3'] != -Fr(mu3['3'][str(nv)].get('(1, 0, 0)', '0')): s_ok = False
    print('(2a) [z4]mu_{4,n} == (n+3)/((2n+3)(4n^2-1)), all n<=', NMAX, ':', p_ok)
    print('(2b) [z2^2]mu_{4,n} == -n/((2n+3)(4n^2-1)), all n<=', NMAX, ':', q_ok)
    print('(2c) [z3]mu_{4,n} == -[z2]mu_{3,n}, all n<=320 :', s_ok)
