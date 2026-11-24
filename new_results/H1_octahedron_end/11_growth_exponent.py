"""
H1: amplitude/exponent fit of the large-order law  a_k/k! ~ Re[ Z k^beta t*^{-k} ],  t* = (-log(n-1) + i pi)/2,
using pairs of consecutive coefficients (assumes the predicted t*; checks that beta is stable).
"""
import json
import mpmath as mp
mp.mp.dps = 60
for n in (2, 3, 4, 5):
    for var in ('eta', 'd'):
        a = [mp.mpf(x) for x in json.load(open(f'hi_coeffs_n{n}_K80.json'))[var]]
        t = (-mp.log(n - 1) + 1j * mp.pi) / 2
        y = [a[k] / mp.factorial(k) for k in range(len(a))]
        # y_k = Re(Z w^k), w = 1/t ; from y_k, y_{k+1}: solve for Z (complex) assuming Z fixed
        w = 1 / t
        Zs = {}
        for k in range(20, 80):
            M = mp.matrix([[mp.re(w ** k), -mp.im(w ** k)], [mp.re(w ** (k + 1)), -mp.im(w ** (k + 1))]])
            sol = mp.lu_solve(M, mp.matrix([y[k], y[k + 1]]))
            Zs[k] = sol[0] + 1j * sol[1]
        betas = []
        for k in (30, 40, 50, 60, 70, 78):
            if k + 1 in Zs and k - 1 in Zs:
                b = (mp.log(abs(Zs[k + 1])) - mp.log(abs(Zs[k - 1]))) / (mp.log(k + 1) - mp.log(k - 1))
                betas.append((k, mp.nstr(b, 6), mp.nstr(abs(Zs[k]) * mp.mpf(k) ** mp.mpf(1.5), 6), mp.nstr(mp.arg(Zs[k]), 6)))
        print(f'n={n} [{var}]  (k, local beta, |Z_k| k^1.5, arg Z_k): {betas}')
