"""
VERIFIER V3 (independent of the solver's explicit beta_k formula and of its 16000-term tail fit).
With f = k1 (first-order coefficient of h_p in eps) and g = F2 = k2 + k1^2/2 (second-order coefficient),
Parseval on [0, pi/2] for cos(4k th) series gives EXACTLY
    sum_{k>=1} (16k^2-1) alpha_k beta_k = 2 avg(f' g' - f g) + alpha_0 beta_0 / 2,   avg = (2/pi) int_0^{pi/2}.
We evaluate this single integral (no series, no tail) and compare with the claimed
    D_3/(4pi^2) = -(15/32) zeta3 + 43 pi^2/384 - (7/48) pi^2 log2 + log2 - 5/16.
Also: the eps^2 analogue sum (16k^2-1) alpha_k^2 = 2 avg(f'^2 - f^2) + alpha_0^2/2 vs 7pi^2/48 - 1,
and the h_p eps-coefficients are re-derived by numerical eps-differentiation of h_p at a point.
"""
import mpmath as mp
mp.mp.dps = 50
pi = mp.pi; lam = mp.log(2); z3 = mp.zeta(3)

def parts(th):
    c, s = mp.cos(th), mp.sin(th)
    C, S = c * c, s * s
    a, b = mp.log(abs(c)), mp.log(abs(s))
    k1 = C * a + S * b
    k2 = C * S * (a - b) ** 2
    F2 = k2 + k1 ** 2 / 2
    s2t, c2t = mp.sin(2 * th), mp.cos(2 * th)
    k1p = s2t * (b - a)                          # derived by hand: k1' = sin2t (log s - log c)
    k2p = s2t * c2t * (a - b) ** 2 - s2t * (a - b)
    F2p = k2p + k1 * k1p
    return k1, k1p, F2, F2p

# check the hand derivatives and the eps-expansion of h_p at a point
th0 = mp.mpf('0.53')
k1, k1p, F2, F2p = parts(th0)
nd1 = mp.diff(lambda t: parts(t)[0], th0); nd2 = mp.diff(lambda t: parts(t)[2], th0)
print("derivative checks:", mp.nstr(k1p - nd1, 3), mp.nstr(F2p - nd2, 3))
hp = lambda e: (abs(mp.cos(th0)) ** (2 / (1 - e)) + abs(mp.sin(th0)) ** (2 / (1 - e))) ** ((1 - e) / 2)
tay = mp.taylor(hp, 0, 2)
print("h_p eps-coefficients vs (k1, F2):", mp.nstr(tay[1] - k1, 3), mp.nstr(tay[2] - F2, 3))

avg = lambda F: 2 / pi * mp.quad(F, [0, pi / 8, pi / 4, 3 * pi / 8, pi / 2])
alpha0 = 4 / pi * mp.quad(lambda t: parts(t)[0], [0, pi / 4, pi / 2])
beta0 = 4 / pi * mp.quad(lambda t: parts(t)[2], [0, pi / 4, pi / 2])
print("alpha_0 =", mp.nstr(alpha0, 30), " (expected 1 - 2 log 2 =", mp.nstr(1 - 2 * lam, 30), ")")
S3 = 2 * avg(lambda t: (lambda q: q[1] * q[3] - q[0] * q[2])(parts(t))) + alpha0 * beta0 / 2
D3c = -mp.mpf(15) / 32 * z3 + 43 * pi ** 2 / 384 - mp.mpf(7) / 48 * pi ** 2 * lam + lam - mp.mpf(5) / 16
print(f"sum (16k^2-1) alpha_k beta_k (integral form) = {mp.nstr(S3, 45)}")
print(f"claimed closed form                           = {mp.nstr(D3c, 45)}   diff {mp.nstr(S3 - D3c, 3)}")
S2 = 2 * avg(lambda t: (lambda q: q[1] ** 2 - q[0] ** 2)(parts(t))) + alpha0 ** 2 / 2
print(f"sum (16k^2-1) alpha_k^2 (integral form) = {mp.nstr(S2, 40)}  vs 7pi^2/48-1 diff {mp.nstr(S2 - (7*pi**2/48 - 1), 3)}")
# PSLQ on S3: is the zeta(3) coefficient -15/32 on a small natural basis?
basis = [S3, z3, pi ** 2, pi ** 2 * lam, lam, 1]
rel = mp.pslq(basis, maxcoeff=10 ** 6, maxsteps=10 ** 6)
print("PSLQ [S3, z3, pi^2, pi^2 log2, log2, 1]:", rel,
      " residual", mp.nstr(abs(mp.fsum(r * b for r, b in zip(rel, basis))), 3) if rel else "")
