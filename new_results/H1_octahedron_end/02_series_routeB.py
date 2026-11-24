"""
H1, octahedron end.  Route B (independent of Route A): cumulant / log-moment expansion.

Representation used here (Sec. 1 of result.md, "centred form"):
  S_n/S_n(1) = (n-1)! eps^{n-1} Gamma(1+a)^n / Gamma(1+(n-1)eps) * E[ f(eta, m) ],
  a = -(n-1) eta / n,  f = sqrt( (1/n) sum_i exp(2 eta m_i) ),  m_i = l_i - mean(l),
  l_i = log G_i, G_i iid Gamma(1+a).   (f is invariant under a common shift of the l_i.)
Expand f in eta as a polynomial in l_1..l_n, then take expectations using the cumulants
  kappa_1 = 0 (shift), kappa_j = psi^{(j-1)}(1+a) = sum_r psi^{(j-1+r)}(1) a^r / r!,
  psi^{(m)}(1) = (-1)^{m+1} m! zeta(m+1).
No finite differences of Gamma and no binomial expansion of sqrt(sum G^{2eta}) are used.
"""
import sys, json, math
from sympy import QQ, Rational, binomial
from sympy.polys.rings import ring

K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
n = int(sys.argv[2]) if len(sys.argv) > 2 else 3

# stage 1: polynomial ring in eta, l1..ln
lnames = [f'l{i}' for i in range(1, n + 1)]
R1, eta1, *ls = ring(','.join(['eta'] + lnames), QQ)


def tr1(p):
    return R1({m: c for m, c in p.items() if m[0] <= K})


mbar = sum(ls) * QQ(1, n)
ms = [li - mbar for li in ls]
y = R1(0)
for mi in ms:
    term = R1(1)
    e = R1(1)
    for j in range(1, K + 1):
        term = tr1(term * (2 * eta1 * mi)) * QQ(1, j)
        e += term
    y += e
y = tr1(y * QQ(1, n) - 1)
# y = O(eta^2)
assert all(m[0] >= 2 for m in y.keys())
f = R1(0)
yk = R1(1)
for k in range(0, K // 2 + 1):
    _b = Rational(binomial(Rational(1, 2), k))
    f += QQ(int(_b.p), int(_b.q)) * yk
    yk = tr1(yk * y)
f = tr1(f)

# stage 2: ring in eta, g, z2..zK
names = ['eta', 'g'] + [f'z{j}' for j in range(2, K + 2)]
R, *gens = ring(','.join(names), QQ)
eta, g = gens[0], gens[1]
z = {j: gens[j] for j in range(2, K + 2)}


def tr(p):
    return R({m: c for m, c in p.items() if m[0] <= K})


def exp_ser(L):
    res = R(1)
    term = R(1)
    for m in range(1, K + 1):
        term = tr(term * L) * QQ(1, m)
        if term == 0:
            break
        res += term
    return res


def logGamma1(c):
    c = QQ(c)
    L = -g * c * eta
    for j in range(2, K + 1):
        L += QQ((-1) ** j, j) * z[j] * c ** j * eta ** j
    return tr(L)


a_c = QQ(-(n - 1), n)  # a = a_c * eta


def psi_der_at_1(mord):
    """psi^{(m)}(1) for m>=1 as ring element"""
    return (-1) ** (mord + 1) * math.factorial(mord) * z[mord + 1]


def kappa(j):
    """kappa_j(eta) = psi^{(j-1)}(1 + a_c eta), j >= 2"""
    s = R(0)
    for r in range(0, K + 1):
        if j - 1 + r + 1 > K + 1:
            break
        s += psi_der_at_1(j - 1 + r) * (a_c * eta) ** r * QQ(1, math.factorial(r))
    return tr(s)


kap = {1: R(0)}
for j in range(2, K + 1):
    kap[j] = kappa(j)
mom = {0: R(1)}
for k in range(1, K + 1):
    s = R(0)
    for j in range(1, k + 1):
        s += math.comb(k - 1, j - 1) * tr(kap[j] * mom[k - j])
    mom[k] = tr(s)

Ef = R(0)
for mono, c in f.items():
    e_eta = mono[0]
    prod = R(QQ(c)) * eta ** e_eta
    for ei in mono[1:]:
        prod = tr(prod * mom[ei])
    Ef += prod
Ef = tr(Ef)

# prefactor (n-1)! eps^{n-1} Gamma(1+a)^n / Gamma(n-(n-1)eta)
pref = exp_ser(tr(n * logGamma1(a_c) - logGamma1(-(n - 1))))
for i in range(1, n):
    geo = R(0)
    for mm in range(K + 1):
        geo += (QQ(n - 1, i) * eta) ** mm
    pref = tr(pref * tr(geo))
for i in range(n - 1):
    pref = tr(pref * (1 - eta))
ser = tr(pref * Ef)

# to d
etad = R(0)
for mm in range(1, K + 1):
    etad += (-1) ** (mm - 1) * eta ** mm
coeffs = {}
for mono, c in ser.items():
    coeffs.setdefault(mono[0], R(0))
    coeffs[mono[0]] += R({(0,) + mono[1:]: c})
sd = R(0)
pw = R(1)
for e in range(K + 1):
    if e > 0:
        pw = tr(pw * etad)
    if e in coeffs:
        sd += tr(coeffs[e] * pw)
sd = tr(sd)
cd = {}
for mono, c in sd.items():
    cd.setdefault(mono[0], R(0))
    cd[mono[0]] += R({(0,) + mono[1:]: c})
ce = coeffs
out = {'eta': [], 'd': []}
print(f'===== Route B, n = {n}, K = {K} =====')
for k in range(K + 1):
    ee = ce.get(k, R(0)).as_expr()
    dd = cd.get(k, R(0)).as_expr()
    print(f'[eta^{k}] = {ee}')
    print(f'[d^{k}]   = {dd}')
    out['eta'].append(str(ee))
    out['d'].append(str(dd))
json.dump(out, open(f'routeB_n{n}_K{K}.json', 'w'), indent=1)
print('saved', f'routeB_n{n}_K{K}.json')
