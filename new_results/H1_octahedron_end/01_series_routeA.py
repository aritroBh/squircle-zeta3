"""
H1, octahedron end.  Route A: exact Taylor coefficients of the cross-polytope-end surface area.

Exact representation (proved in result.md, Sec. 1):
  S_n(p) = 2^n eps^{n-1} Gamma(eps)^n / Gamma(1+(n-1)eps) * E[ sqrt( sum_i G_i^{2 eta} ) ],
  G_i iid Gamma(eps,1),  eps = 1/p,  eta = 1 - eps.
Write sqrt(sum X_i) = sqrt(n) * sqrt(1+x), x = (1/n) sum (X_i - 1), expand sqrt(1+x) binomially;
E[(G^{2eta}-1)^k] = Phi_k / Gamma(eps),  Phi_k = sum_j C(k,j)(-1)^{k-j} Gamma(1+(2j-1)eta).
Truncating at k <= K gives the exact eta-expansion through eta^K (remainder O(eta^{K+1})).
  S_n/S_n(1) = (1-eta)^{n-1} (n-1)!/Gamma(n-(n-1)eta) * sum_k binom(1/2,k) n^{-k} T_k,
  T_k = sum over compositions (k_1..k_n) of k: multinomial * prod Phi_{k_i}.
Gamma(1+x) = exp(-g x + sum_{j>=2} (-1)^j z_j x^j / j), g = Euler gamma, z_j = zeta(j) (symbols).

Output: coefficients c_k of S_n(1+d)/S_n(1) = sum c_k d^k as exact polynomials in g, z_j.
"""
import sys, json, itertools, math
from sympy import Rational, binomial, factorial, QQ, symbols, Symbol, sympify, expand, nsimplify
from sympy.polys.rings import ring

K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
NS = [int(a) for a in sys.argv[2:]] if len(sys.argv) > 2 else [2, 3, 4]

names = ['eta', 'g'] + [f'z{j}' for j in range(2, K + 1)]
R, *gens = ring(','.join(names), QQ)
eta, g = gens[0], gens[1]
z = {j: gens[j] for j in range(2, K + 1)}  # gens[j] is z_j because gens[2] = z2


def trunc(p, K=K):
    """drop terms of eta-degree > K"""
    return R({m: c for m, c in p.items() if m[0] <= K})


def mul(a, b):
    return trunc(a * b)


def exp_ser(L):
    """exp(L) truncated, L has no constant term"""
    res = R(1)
    term = R(1)
    for m in range(1, K + 1):
        term = mul(term, L) * QQ(1, m)
        if term == 0:
            break
        res += term
    return res


def logGamma1(c):
    """log Gamma(1 + c*eta) as series in eta, c rational"""
    c = QQ(c)
    L = -g * c * eta
    for j in range(2, K + 1):
        L += QQ((-1) ** j, j) * z[j] * c ** j * eta ** j
    return trunc(L)


def Gamma1(c):
    return exp_ser(logGamma1(c))


def invGamma1(c):
    return exp_ser(-logGamma1(c))


def compositions(k, n):
    if n == 1:
        yield (k,)
        return
    for i in range(k + 1):
        for rest in compositions(k - i, n - 1):
            yield (i,) + rest


def series_in_eta(n):
    Phi = []
    for k in range(K + 1):
        s = R(0)
        for j in range(k + 1):
            s += math.comb(k, j) * (-1) ** (k - j) * Gamma1(2 * j - 1)
        Phi.append(trunc(s))
    # sanity: Phi_k = O(eta^k)
    for k, ph in enumerate(Phi):
        assert all(m[0] >= k for m in ph.keys()), k
    E = R(0)
    for k in range(K + 1):
        Tk = R(0)
        for comp in compositions(k, n):
            mult = math.factorial(k)
            for ki in comp:
                mult //= math.factorial(ki)
            prod = R(1)
            for ki in comp:
                prod = mul(prod, Phi[ki])
            Tk += mult * prod
        _b = Rational(binomial(Rational(1, 2), k)); bk = QQ(int(_b.p), int(_b.q))
        E += bk * QQ(1, n ** k) * Tk
    E = trunc(E)
    # prefactor (1-eta)^{n-1} (n-1)!/Gamma(n-(n-1)eta)
    #   Gamma(n - (n-1)eta) = Gamma(1-(n-1)eta) * prod_{i=1}^{n-1} (i - (n-1) eta)
    #   (n-1)!/that = invGamma1(-(n-1)) * prod 1/(1 - (n-1)eta/i)
    pref = invGamma1(-(n - 1))
    for i in range(1, n):
        geo = R(0)
        for m in range(K + 1):
            geo += (QQ(n - 1, i) * eta) ** m
        pref = mul(pref, trunc(geo))
    for i in range(n - 1):
        pref = mul(pref, 1 - eta)
    # 1 - eta = eps; also Gamma(eps)^n cancels exactly (see docstring)
    # NOTE: the Gamma(eps)^n / Gamma(eps)^n cancellation: prefactor has Gamma(eps)^n, E[...] has 1/Gamma(eps)^n
    # per factor, but Phi_0 = Gamma(1-eta) = Gamma(eps) carries the Gamma(eps) back.  So:
    #   S_n/S_n(1) = (n-1)! eps^{n-1}/Gamma(1+(n-1)eps) * sum_k binom n^{-k} T_k  with T_k built from Phi's.
    return trunc(mul(pref, E))


def to_d(ser_eta):
    """substitute eta = d/(1+d) and re-expand in d (reuse the variable eta as d)"""
    # eta(d) = d - d^2 + d^3 - ...
    etad = R(0)
    for m in range(1, K + 1):
        etad += (-1) ** (m - 1) * eta ** m
    out = R(0)
    powk = R(1)
    # collect coefficient polynomials in eta
    coeffs = {}
    for mono, c in ser_eta.items():
        e = mono[0]
        rest = (0,) + mono[1:]
        coeffs.setdefault(e, R(0))
        coeffs[e] += R({rest: c})
    for e in range(K + 1):
        if e > 0:
            powk = mul(powk, etad)
        if e in coeffs:
            out += mul(coeffs[e], powk)
    return trunc(out)


def coeff_list(ser):
    cl = [R(0)] * (K + 1)
    for mono, c in ser.items():
        e = mono[0]
        cl[e] += R({(0,) + mono[1:]: c})
    return cl


if __name__ == '__main__':
    out = {}
    zsyms = symbols(names)
    for n in NS:
        se = series_in_eta(n)
        sd = to_d(se)
        ce = coeff_list(se)
        cd = coeff_list(sd)
        print(f'===== n = {n}  (S_n(1) = 2^n sqrt(n)/(n-1)!) =====')
        out[n] = {'eta': [], 'd': []}
        for k in range(K + 1):
            ex_e = ce[k].as_expr()
            ex_d = cd[k].as_expr()
            gfree = 'g' not in str(ex_d)
            print(f'[eta^{k}] = {ex_e}')
            print(f'[d^{k}]   = {ex_d}    (gamma-free: {gfree})')
            out[n]['eta'].append(str(ex_e))
            out[n]['d'].append(str(ex_d))
    fname = f'routeA_n{"-".join(map(str, NS))}_K{K}.json'
    with open(fname, 'w') as f:
        json.dump(out, f, indent=1)
    print('saved', fname)
