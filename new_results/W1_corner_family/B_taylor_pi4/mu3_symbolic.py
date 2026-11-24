"""mu3_symbolic.py -- PROVED (exact symbolic algebra over Q(n)).
Closed formula for the flat-angle Taylor coefficients mu_{k,n} (k <= 3) of the corner coefficients,
   M_k(alpha) = sum_{n>=1} mu_{k,n} c^{2n},  c = cos(alpha/2),
as an exact polynomial in
   h1,h2,h3 = H_{n-1}, H^(2)_{n-1}, H^(3)_{n-1},   g1,g2,g3 = H_{2n-1}, H^(2)_{2n-1}, H^(3)_{2n-1},  z2, z3
with coefficients in Q(n), valid for every integer n >= 1.

Derivation (G1 Theorem A; re-implemented independently):
  mu_{k,n} = -(1/(n(2n-1))) sum_{(j,i,l)} e^{(k)}_{j,i,l}(n) (n)_j/(2n)_j D_{i,l}(n+j, n),
  [u^k] E_n = sum e_{j,i,l}(n) y^j a^i b^l,  E_n = exp(uA)(1+y(e^{ul}-1))^{1-2n},
  A = -(n-1)a + n b, l = b - a, a = log y, b = log(1-y);
  D_{i,l}(p,q) = d_p^i d_q^l B(p,q)/B(p,q): moments of the cumulants
  kappa_{r,0} = psi^(r-1)(p) - psi^(r-1)(p+q), kappa_{0,m} = psi^(m-1)(q) - psi^(m-1)(p+q),
  kappa_{r,m} = -psi^(r+m-1)(p+q) (r,m >= 1); moment recursion
  D_{I,L} = sum_{r=1}^{I} sum_{m=0}^{L} C(I-1,r-1) C(L,m) kappa_{r,m} D_{I-r,L-m},
  D_{0,L} = sum_{m=1}^{L} C(L-1,m-1) kappa_{0,m} D_{0,L-m}.
  psi(N)+gamma = H_{N-1}, psi'(N) = z2 - H^(2)_{N-1}, psi''(N) = -2 z3 + 2 H^(3)_{N-1}.
Writes mu_formula.pkl and the Z-sum canonical form mu_zsum.pkl (input to b3_allN_proof.py).
"""
import sys, json, pickle
from fractions import Fraction as Fr
from math import comb
import sympy as sp
from sympy import QQ
from sympy.polys.rings import ring
sys.path.insert(0, '.')
import zsum
from zsum import K, nn

KMAX = 3
NCHK = int(sys.argv[1]) if len(sys.argv) > 1 else 320
u, y, a, b, n = sp.symbols('u y a b n')

A = -(n - 1)*a + n*b
l = b - a
logE = u*A + (1 - 2*n)*sp.log(1 + y*(sp.exp(u*l) - 1))
E_ser = sp.series(sp.exp(sp.series(logE, u, 0, KMAX + 1).removeO()), u, 0, KMAX + 1).removeO()
EK = {k: sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b) for k in range(KMAX + 1)}

R, h1, h2, h3, g1, g2, g3, z2, z3 = ring('h1,h2,h3,g1,g2,g3,z2,z3', K)
NAMES = ['h1', 'h2', 'h3', 'g1', 'g2', 'g3', 'z2', 'z3']


def toK(expr):
    P = sp.Poly(sp.expand(expr), n)
    out = K(0)
    for (e,), c in P.terms():
        out += K(QQ(int(sp.Rational(c).p), int(sp.Rational(c).q)))*nn**e
    return out


def psiR(m, which, j):
    H = {1: h1, 2: h2, 3: h3} if which == 'n' else {1: g1, 2: g2, 3: g3}
    base = nn if which == 'n' else 2*nn
    r = m + 1
    harm = H[r] + sum((K(1)/(base + i)**r for i in range(j)), K(0))
    if m == 0:
        return harm
    if m == 1:
        return z2 - harm
    if m == 2:
        return -2*z3 + 2*harm
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


def eval_formula(F, nv):
    vals = [Hf(1, nv - 1), Hf(2, nv - 1), Hf(3, nv - 1), Hf(1, 2*nv - 1), Hf(2, 2*nv - 1), Hf(3, 2*nv - 1)]
    out = {'1': Fr(0), 'z2': Fr(0), 'z3': Fr(0)}
    for mon, cf in F.terms():
        v = evalK(cf, nv)
        for x, e in zip(vals, mon[:6]):
            v *= x**e
        key = 'z2' if mon[6] else ('z3' if mon[7] else '1')
        assert mon[6] + mon[7] <= 1
        out[key] += v
    return out


def mu_exact(k, nv):
    """independent exact evaluation with numeric polygamma values at integer arguments."""
    def P(m, N):
        if m == 0:
            return (Hf(1, N - 1), Fr(0), Fr(0))
        if m == 1:
            return (-Hf(2, N - 1), Fr(1), Fr(0))
        return (2*Hf(3, N - 1), Fr(0), Fr(-2))
    def add(x, y_, sg=1):
        return tuple(x[i] + sg*y_[i] for i in range(3))
    def mul(x, y_):
        assert x[1]*y_[1] == 0 and x[1]*y_[2] == 0 and x[2]*y_[1] == 0 and x[2]*y_[2] == 0
        return (x[0]*y_[0], x[0]*y_[1] + x[1]*y_[0], x[0]*y_[2] + x[2]*y_[0])
    tot = (Fr(0), Fr(0), Fr(0))
    for (j, i, ll), co in EK[k].terms():
        cv = sp.Rational(co.subs(n, nv)); cv = Fr(int(cv.p), int(cv.q))
        p, q = nv + j, nv
        kap = {}
        for r in range(0, KMAX + 1):
            for m in range(0, KMAX + 1 - r):
                if r + m == 0:
                    continue
                if m == 0:
                    kap[(r, m)] = add(P(r - 1, p), P(r - 1, p + q), -1)
                elif r == 0:
                    kap[(r, m)] = add(P(m - 1, q), P(m - 1, p + q), -1)
                else:
                    x = P(r + m - 1, p + q); kap[(r, m)] = (-x[0], -x[1], -x[2])
        memo = {}
        def Dm(I, Lq):
            if (I, Lq) in memo:
                return memo[(I, Lq)]
            if I == 0 and Lq == 0:
                res = (Fr(1), Fr(0), Fr(0))
            elif I > 0:
                res = (Fr(0), Fr(0), Fr(0))
                for r in range(1, I + 1):
                    for m in range(0, Lq + 1):
                        res = add(res, mul((Fr(comb(I - 1, r - 1)*comb(Lq, m)), Fr(0), Fr(0)), mul(kap[(r, m)], Dm(I - r, Lq - m))))
            else:
                res = (Fr(0), Fr(0), Fr(0))
                for m in range(1, Lq + 1):
                    res = add(res, mul((Fr(comb(Lq - 1, m - 1)), Fr(0), Fr(0)), mul(kap[(0, m)], Dm(0, Lq - m))))
            memo[(I, Lq)] = res
            return res
        pr = Fr(1)
        for tt in range(j):
            pr *= Fr(nv + tt, 2*nv + tt)
        tot = add(tot, mul((cv*pr, Fr(0), Fr(0)), Dm(i, ll)))
    f = Fr(-1, nv*(2*nv - 1))
    return {'1': f*tot[0], 'z2': f*tot[1], 'z3': f*tot[2]}


def to_zsum(F):
    """map the polynomial in h,g (coefficients Q(n)) to Z-sum canonical form, separately per zeta part."""
    sub = [zsum.harm_n_minus_1(1), zsum.harm_n_minus_1(2), zsum.harm_n_minus_1(3),
           zsum.harm_2n_minus_1(1), zsum.harm_2n_minus_1(2), zsum.harm_2n_minus_1(3)]
    parts = {'1': {}, 'z2': {}, 'z3': {}}
    for mon, cf in F.terms():
        E = {(): K(1)}
        for X, e in zip(sub, mon[:6]):
            for _ in range(e):
                E = zsum.mul(E, X)
        key = 'z2' if mon[6] else ('z3' if mon[7] else '1')
        for w, c in E.items():
            zsum.add_into(parts[key], w, cf*c)
    return parts


if __name__ == '__main__':
    res = {}
    for k in range(KMAX + 1):
        F = mu_formula(k)
        res[k] = F
        print(f'mu_{{{k},n}}  ({len(F.terms())} monomials):')
        for mon, cf in sorted(F.terms(), key=lambda x: (sum(x[0]), x[0])):
            name = '*'.join(f'{v}^{e}' if e > 1 else v for v, e in zip(NAMES, mon) if e) or '1'
            print(f'   [{name}] : {sp.factor(K.to_sympy(cf) if hasattr(K, "to_sympy") else cf.as_expr())}')
        sys.stdout.flush()
    # exact checks
    tab = json.load(open('research/10_frontier/G1_zeta3_any_angle/flat_exact_N320_K3.json'))
    keymap = {'1': '(0, 0, 0)', 'z2': '(1, 0, 0)', 'z3': '(0, 1, 0)'}
    for k in range(KMAX + 1):
        ok_ind = True; ok_g1 = True
        for nv in range(1, NCHK + 1):
            fv = eval_formula(res[k], nv)
            if nv <= 40:
                if fv != mu_exact(k, nv):
                    ok_ind = False; print('MISMATCH formula vs independent exact', k, nv)
            g = tab[str(k)][str(nv)]
            if any(kk not in keymap.values() for kk in g):
                ok_g1 = False
            if {kk: Fr(g.get(keymap[kk], '0')) for kk in keymap} != fv:
                ok_g1 = False; print('MISMATCH formula vs G1 table', k, nv)
        print(f'k={k}: formula == independent exact evaluation (n<=40): {ok_ind};  == G1 table (n<={NCHK}): {ok_g1}')
    print('mu_{3,1} =', mu_exact(3, 1))
    print('mu_{3,2} =', mu_exact(3, 2))
    # zeta parts of mu_{3,n} in closed form
    F3 = res[3]
    zc = {'z2': R(0), 'z3': R(0)}
    for mon, cf in F3.terms():
        if mon[6]:
            zc['z2'] += R({tuple(list(mon[:6]) + [0, 0]): cf})
        if mon[7]:
            zc['z3'] += R({tuple(list(mon[:6]) + [0, 0]): cf})
    print('zeta(3)-part of mu_{3,n}:', zc['z3'])
    print('zeta(2)-part of mu_{3,n}:', zc['z2'])
    pickle.dump({k: [(mon, str(cf.as_expr())) for mon, cf in res[k].terms()] for k in res}, open('mu_formula.pkl', 'wb'))
    Z = to_zsum(F3)
    pickle.dump({part: [(w, str(c.as_expr())) for w, c in E.items()] for part, E in Z.items()}, open('mu_zsum.pkl', 'wb'))
    print('Z-sum canonical form of mu_{3,n}: words per part =', {p: len(E) for p, E in Z.items()})
    # check the Z-sum form numerically-exactly
    ok = True
    for nv in range(2, 30):
        for part, E in Z.items():
            v = sum((evalK(c, nv)*zsum.Z_val(w, 2*nv) for w, c in E.items()), Fr(0))
            if v != eval_formula(F3, nv)[part]:
                ok = False; print('zsum mismatch', part, nv)
    print('Z-sum form == formula for 2<=n<30:', ok)
