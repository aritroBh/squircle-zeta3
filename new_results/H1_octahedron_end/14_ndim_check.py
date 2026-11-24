"""
H1: numerical check of the n-dimensional expansion (n = 3, 4, 5) without series expansion.
  S_n(p) = 2^n eps^{n-1} Gamma(eps)^n / Gamma(1+(n-1)eps) * E[ sqrt(sum_i G_i^{2 eta}) ],  G_i iid Gamma(eps)
  sqrt(X) = (1/(2 sqrt(pi))) Int_0^inf (1 - e^{-tX}) t^{-3/2} dt  =>
  E[sqrt(sum G_i^{2eta})] = (1/(2 sqrt(pi))) Int_0^inf (1 - h(t)^n) t^{-3/2} dt,
  h(t) = E[exp(-t G^{2eta})] = (1/Gamma(eps)) Int_0^inf exp(-t g^{2eta}) g^{eps-1} e^{-g} dg.
Compared with the truncated Taylor series from hi_coeffs_n{n}_K80.json (orders <= 12 exact-checked).
For n = 3 also compared with the direct graph quadrature (direct_S_dps40.json).
"""
import json
import mpmath as mp

mp.mp.dps = 32


def S_gamma(n, p):
    p = mp.mpf(p)
    eps = 1 / p
    eta = 1 - eps
    Ge = mp.gamma(eps)

    def h(t):
        return mp.quad(lambda g: mp.exp(-t * g ** (2 * eta) - g) * g ** (eps - 1), [0, 1, mp.inf]) / Ge

    E = mp.quad(lambda t: (1 - h(t) ** n) * t ** mp.mpf(-1.5), [0, 1, mp.inf]) / (2 * mp.sqrt(mp.pi))
    return 2 ** n * eps ** (n - 1) * Ge ** n / mp.gamma(1 + (n - 1) * eps) * E


direct = json.load(open('direct_S_dps40.json'))
for n in (3, 4, 5):
    c = [mp.mpf(x) for x in json.load(open(f'hi_coeffs_n{n}_K80.json'))['d']]
    S1 = 2 ** n * mp.sqrt(n) / mp.factorial(n - 1)
    for dstr in ('-0.01', '0.01', '-0.001', '0.001'):
        d = mp.mpf(dstr)
        val = S_gamma(n, 1 + d)
        ser12 = S1 * sum(c[k] * d ** k for k in range(13))
        ser3 = S1 * sum(c[k] * d ** k for k in range(4))
        extra = ''
        if n == 3:
            extra = f'  |gamma-rep - direct graph quadrature| = {mp.nstr(abs(val - mp.mpf(direct[dstr]["N1"])), 3)}'
        print(f'n={n} d={dstr:>6}: S = {mp.nstr(val, 28)}  |S - series_12| = {mp.nstr(abs(val - ser12), 3)}  '
              f'(S - series_3)/(S1 d^4) = {mp.nstr((val - ser3) / (S1 * d ** 4), 8)} vs c4 = {mp.nstr(c[4], 8)}' + extra, flush=True)
