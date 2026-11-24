"""
H1: final numbers and the general-n statements.
 (a) log B_n at eta^4:  -(n-1)/n^2 zeta(2)^2 + (n-1)(n^2-3n+4)/n^3 zeta(4)   (check against exact data, n = 2..12)
 (b) explicit d-coefficients for general n:
        c1 = sigma1 = (n-1)(H_{n-1}-1),  c2 = (sigma1^2+sigma2)/2 + A2,  c3 = sigma1^3/6 + sigma1(sigma2/2+A2) + sigma3/3 - 2A2 + A3,
        sigma_k = sum_{i=1}^{n-2} ((n-1-i)/i)^k,  A2 = -(n-1)(n^2-2n-1)/(2n) zeta2,  A3 = -(n^2-1)((n-1)^3-2)/(3n^2) zeta3
     (check against exact Route A d-coefficients for every available n)
 (c) n = 3: exact S_k and 60-digit values, k <= 8.
"""
import json, glob
import sympy as sp
import mpmath as mp

z = {j: sp.Symbol(f'z{j}') for j in range(2, 30)}
even = {z[2 * j]: sp.zeta(2 * j) for j in range(1, 14)}
LB = json.load(open('logB_coeffs.json'))
print('(a) eta^4 coefficient of log B_n')
for n, rows in sorted(LB.items(), key=lambda t: int(t[0])):
    n = int(n)
    got = sp.sympify(rows[2])  # rows start at eta^2
    pred = -sp.Rational(n - 1, n * n) * z[2] ** 2 + sp.Rational((n - 1) * (n * n - 3 * n + 4), n ** 3) * z[4]
    print(f'  n={n:>2}: got - pred = {sp.simplify(sp.expand(got - pred).subs(even))}')

print('\n(b) general-n d-coefficients c1, c2, c3')
data = {}
for f in sorted(glob.glob('routeA_n*_K*.json')):
    for n, v in json.load(open(f)).items():
        if int(n) not in data or len(data[int(n)]) < len(v['d']):
            data[int(n)] = v['d']
for n in sorted(data):
    sig = lambda k: sum(sp.Rational(n - 1 - i, i) ** k for i in range(1, n - 1))
    H = sum(sp.Rational(1, i) for i in range(1, n))
    A2 = -sp.Rational((n - 1) * (n * n - 2 * n - 1), 2 * n) * z[2]
    A3 = -sp.Rational((n * n - 1) * ((n - 1) ** 3 - 2), 3 * n * n) * z[3]
    c1 = sig(1)
    c2 = (sig(1) ** 2 + sig(2)) / 2 + A2
    c3 = sig(1) ** 3 / 6 + sig(1) * (sig(2) / 2 + A2) + sig(3) / 3 - 2 * A2 + A3
    d = [sp.sympify(s) for s in data[n]]
    chk = [sp.simplify(sp.expand(x - y).subs(even)) for x, y in ((d[1], c1), (d[2], c2), (d[3], c3))]
    print(f'  n={n:>2}: c1 = {c1} (=(n-1)(H_(n-1)-1): {sp.simplify(c1 - (n - 1) * (H - 1)) == 0}), '
          f'c2 = {sp.expand(c2.subs(even))}, c3 = {sp.expand(c3.subs(even))};  exact-data checks {chk}')

print('\n(c) n = 3 (octahedron): S(1+d) = 4 sqrt3 * sum c_k d^k')
mp.mp.dps = 70
d3 = data[3]
for k in range(0, 9):
    ek = sp.expand(sp.sympify(d3[k]).subs(even))
    val = mp.mpf(str(sp.N(ek.subs({z[j]: sp.zeta(j) for j in range(2, 30)}), 75)))
    print(f'  c_{k} = {ek}')
    print(f'       = {mp.nstr(val, 60)};   S_{k} = 4 sqrt3 c_{k} = {mp.nstr(4 * mp.sqrt(3) * val, 60)}')
