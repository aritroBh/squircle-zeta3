"""PIPDIV script 3: numerical checks of the identities used in the proof.
 (a) Laplace form:  E/4 = h sum_{n>=1} binom(h,n) int_0^oo e^{-n y/h} e^{(n+1)y} (1+e^y)^{-1-h} dy
     against the s-integral and against pi_10 = 3.73396080922282829958... (D12 table).
 (b) branch check: A(y,h) = e^{2y}(1+e^y)^{-1-h} = y0^{-1-h}(1-y/y0)^{-1-h} Phi(y-y0,h) near y = 0,
     y0 = +-i pi, complex h, Phi(u,h) = e^{2u} ((e^u-1)/u)^{-1-h} (log of (e^u-1)/u taken analytic, = principal here).
 (c) |Gamma(M+h)/Gamma(M)| <= 1.4 M^{|h|} for M >= 1, |h| <= 1/4 (grid).
 (d) Lemma 1 at order h^{K}: (E/4 - sum_{k<K} T-coefficients) / h^{K+1} -> N_K, using the pd1 table."""
import sys
import mpmath as mp
mp.mp.dps = 30

def E4_s(h):
    k = h / (1 - h)
    f = lambda s: h * k * mp.exp(-h * k * s) * (1 + mp.exp(-k * s)) ** (-1 - h) * ((1 + mp.exp(-s)) ** h - 1)
    return mp.quad(f, [0, 1, 10, 50, mp.inf])
def E4_laplace(h):
    def term(n):
        n = int(n)
        g = lambda y: mp.exp(-n * y / h + (n + 1) * y) * (1 + mp.exp(y)) ** (-1 - h)
        return h * mp.binomial(h, n) * mp.quad(g, [0, h / n, 10 * h / n, mp.inf])
    return mp.nsum(term, [1, mp.inf])
for p in (10, 4):
    h = mp.mpf(1) / p
    a, b = E4_s(h), E4_laplace(h)
    print('(a) p=%d: s-form %s  Laplace form %s  diff %s' % (p, mp.nstr(a, 20), mp.nstr(b, 20), mp.nstr(a - b, 3)))
    if p == 10:
        print('    4(2^-h + E/4) = %s   (pi_10 = 3.73396080922282829958)' % mp.nstr(4 * (2 ** -h + a), 21))

phi_log = lambda u: mp.log((mp.exp(u) - 1) / u)
for y0 in (mp.mpc(0, mp.pi), mp.mpc(0, -mp.pi)):
    for h in (mp.mpc(0.2, 0.1), mp.mpc(-0.15, -0.2)):
        for y in (mp.mpf(0), mp.mpc(0.3, 0.4), mp.mpc(-0.5, 0.2)):
            u = y - y0
            lhs = mp.exp(2 * y) * mp.exp(-(1 + h) * mp.log(1 + mp.exp(y)))
            rhs = mp.exp(-(1 + h) * mp.log(y0)) * mp.exp(-(1 + h) * mp.log(1 - y / y0)) * mp.exp(2 * u - (1 + h) * phi_log(u))
            print('(b) y0=%s h=%s y=%s  |A - (y0-y)^{-1-h}Phi| = %s' % (mp.nstr(y0, 3), mp.nstr(h, 3), mp.nstr(y, 3), mp.nstr(abs(lhs - rhs), 3)))

worst = 0
for M in [1 + i / 20 for i in range(0, 80)] + [5, 10, 100, 1000, 10**6]:
    for j in range(64):
        h = 0.25 * mp.expjpi(mp.mpf(j) / 32)
        worst = max(worst, abs(mp.gamma(M + h) / mp.gamma(M)) / mp.mpf(M) ** 0.25)
print('(c) max |Gamma(M+h)/Gamma(M)| / M^{1/4} on grid = %s (claimed <= 1.4)' % mp.nstr(worst, 6))

if len(sys.argv) > 1:
    mp.mp.dps = 40
    N = [mp.mpf(l) for l in open(sys.argv[1])][:8]
    for p in (200, 400, 800):
        h = mp.mpf(1) / p
        r = 2 ** -h + E4_s(h) - 1 - sum(N[k] * h ** (k + 1) for k in range(5))
        print('(d) p=%d: (pi_p/4 - 1 - sum_{k<5} N_k h^{k+1})/h^6 = %s   N_5 = %s' % (p, mp.nstr(r / h ** 6, 10), mp.nstr(N[5], 10)))
