"""
VERIFIER V2 (independent of the solver's support-function / Beta-derivative route).
Circle end, eps = 2/p - 1.  P(eps) computed by the POLAR ARCLENGTH of the radial function
    rho = (c^p + s^p)^(-1/p),  P = 8 int_0^{pi/4} sqrt(rho^2 + rho'^2) dth,
at Chebyshev nodes in eps in [-r, r]; Taylor coefficients P_0..P_5 from the Chebyshev interpolant.
Then D = P^2 - 4 pi A with the exact area A = 4 Gamma((3+eps)/2)^2/Gamma(2+eps).
Compared with the solver's claimed closed forms for P_0..P_3, D_2/(2pi^2) and D_3/(4pi^2).
"""
import mpmath as mp
mp.mp.dps = 60
pi = mp.pi; lam = mp.log(2); z3 = mp.zeta(3)

def P_polar(p):
    def f(th):
        c, s = mp.cos(th), mp.sin(th)
        S = c ** p + s ** p
        rho = S ** (-1 / p)
        drho = (c ** (p - 1) * s - s ** (p - 1) * c) * S ** (-1 / p - 1)
        return mp.sqrt(rho * rho + drho * drho)
    return 8 * mp.quad(f, [0, pi / 16, pi / 8, pi / 4])

print("sanity P(2) - 2pi =", mp.nstr(P_polar(mp.mpf(2)) - 2 * pi, 3))
N = 44; r = mp.mpf('0.1')
nodes = [r * mp.cos(pi * (j + mp.mpf(1) / 2) / N) for j in range(N)]
vals = [P_polar(2 / (1 + e)) for e in nodes]
# Chebyshev coefficients
cheb = []
for m in range(N):
    cm = 2 * mp.fsum(vals[j] * mp.cos(m * pi * (j + mp.mpf(1) / 2) / N) for j in range(N)) / N
    cheb.append(cm)
cheb[0] /= 2
print("last Chebyshev coefficients (decay check):", [mp.nstr(abs(cheb[m]), 3) for m in range(N - 6, N)])
# convert sum c_m T_m(x), x = eps/r  to a power series in eps
# T_m coefficients via recurrence (as polynomials in x)
T = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]]
for m in range(2, N):
    a = [mp.mpf(0)] + [2 * c for c in T[m - 1]]
    b = T[m - 2] + [mp.mpf(0)] * (len(a) - len(T[m - 2]))
    T.append([x - y for x, y in zip(a, b)])
Pk = []
for n in range(6):
    Pk.append(mp.fsum(cheb[m] * T[m][n] for m in range(n, N) if n < len(T[m])) / r ** n)
claim = {
    0: 2 * pi,
    1: pi * (1 - 2 * lam),
    2: pi * (lam ** 2 - lam - mp.mpf(3) / 4 + 11 * pi ** 2 / 96),
    3: pi * (-mp.mpf(23) / 32 * z3 - mp.mpf(11) / 96 * pi ** 2 * lam - lam ** 3 / 3 + lam ** 2 / 2 + mp.mpf(3) / 4 * lam
             + mp.mpf(1) / 16 + 37 * pi ** 2 / 384),
}
for n in range(6):
    s_ = f"P_{n} = {mp.nstr(Pk[n], 40)}"
    if n in claim:
        s_ += f"   claim diff {mp.nstr(Pk[n] - claim[n], 3)}"
    print(s_)
# area series (exact Gamma) and deficit
Aser = mp.taylor(lambda e: 4 * mp.gamma((3 + e) / 2) ** 2 / mp.gamma(2 + e), 0, 5)
Dser = [mp.fsum(Pk[i] * Pk[n - i] for i in range(n + 1)) - 4 * pi * Aser[n] for n in range(6)]
print("D_0, D_1 =", mp.nstr(Dser[0], 3), mp.nstr(Dser[1], 3))
D2c = 7 * pi ** 2 / 48 - 1
D3c = -mp.mpf(15) / 32 * z3 + 43 * pi ** 2 / 384 - mp.mpf(7) / 48 * pi ** 2 * lam + lam - mp.mpf(5) / 16
print(f"D_2/(2pi^2) = {mp.nstr(Dser[2] / (2 * pi ** 2), 40)}  claim 7pi^2/48-1 diff {mp.nstr(Dser[2]/(2*pi**2) - D2c, 3)}")
print(f"D_3/(4pi^2) = {mp.nstr(Dser[3] / (4 * pi ** 2), 40)}  claim diff {mp.nstr(Dser[3]/(4*pi**2) - D3c, 3)}")
# compare also with the stored circle_end exact_values (P_k there)
import os
fn = 'research/05_final_verify/circle_end/exact_values.txt'
if os.path.exists(fn):
    for line in open(fn):
        w = line.split()
        if len(w) == 3 and w[0] == 'P' and int(w[1]) < 6:
            print(f"   circle_end stored P_{w[1]} - this run = {mp.nstr(mp.mpf(w[2]) - Pk[int(w[1])], 3)}")
