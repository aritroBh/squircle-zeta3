"""
F1 step 6: (a) distributional form of the diamond theorem for a non-symmetric smooth test function phi:
      int_0^{2pi} (h_p - h_1) phi dtheta = (sqrt2/4)(zeta(2) tau^2 + zeta(3) tau^3) sum_{j=0}^{3} phi(pi/4 + j pi/2) + O(tau^4)
(b) figure: the boundary layer and the normalized Fourier coefficients.
"""
import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

mp.mp.dps = 40
pi = mp.pi; s2 = mp.sqrt(2)
phi = lambda t: mp.exp(mp.sin(t) + mp.cos(3 * t) / 2)   # no symmetry at all

def h(tau, t):
    q = 1 / tau
    return (abs(mp.cos(t)) ** q + abs(mp.sin(t)) ** q) ** tau
def h1(t):
    return max(abs(mp.cos(t)), abs(mp.sin(t)))

def lhs(tau):
    tot = 0
    for j in range(4):
        tj = pi / 4 + j * pi / 2
        pts = [tj - pi / 4, tj - 8 * tau, tj - 2 * tau, tj, tj + 2 * tau, tj + 8 * tau, tj + pi / 4]
        tot += mp.quad(lambda t: (h(tau, t) - h1(t)) * phi(t), pts)
    return tot
comb = mp.fsum(phi(pi / 4 + j * pi / 2) for j in range(4))
print("tau, [LHS - (sqrt2/4)(z2 tau^2 + z3 tau^3) comb]/tau^4")
for tau in [mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]:
    r = (lhs(tau) - s2 / 4 * (mp.zeta(2) * tau ** 2 + mp.zeta(3) * tau ** 3) * comb) / tau ** 4
    print(mp.nstr(tau, 3), mp.nstr(r, 12))
print("(bounded => the tau^2 and tau^3 terms are exactly the Dirac comb at the edge normals, no phi' terms)")

# ---------- figure ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
th = np.linspace(0, np.pi / 2, 2001)
for tau in [0.15, 0.08, 0.04]:
    q = 1 / tau
    hp = (np.abs(np.cos(th)) ** q + np.abs(np.sin(th)) ** q) ** tau
    ax[0].plot(th, (hp - np.maximum(np.abs(np.cos(th)), np.abs(np.sin(th)))) / tau ** 2, label=f"tau = {tau}")
ax[0].axvline(np.pi / 4, color='gray', lw=0.6)
ax[0].set_xlabel("theta"); ax[0].set_ylabel("(h_p - h_1)/tau^2")
ax[0].set_title("support function: rounding concentrates at the edge normal")
ax[0].legend()

def a_minus_A(k, tau):
    Gk = lambda x: (lambda th_: mp.cos(4 * k * th_) * mp.cos(th_) ** 2 * mp.sin(th_))(mp.atan(mp.e ** (-x)))
    f = lambda s: (mp.e ** (tau * mp.log1p(mp.e ** (-s))) - 1) * Gk(tau * s)
    return 8 / pi * tau * mp.quad(f, [0, 2, 8, 30, 120, mp.inf])
mp.mp.dps = 20
ks = list(range(0, 41))
for tau in [mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.02')]:
    vals = [float((-1) ** k * pi * a_minus_A(k, tau) / (s2 * tau ** 2)) for k in ks]
    ax[1].plot(ks, vals, 'o-', ms=3, label=f"tau = {tau}")
    kk = np.array(ks[1:], dtype=float); t = float(tau)
    lead = (1 / (4 * kk ** 2) - np.pi * t / (2 * kk * np.sinh(2 * np.pi * kk * t))) / t ** 2
    ax[1].plot(kk, lead, 'k:', lw=0.8)
ax[1].axhline(float(mp.zeta(2)), color='gray', lw=0.6)
ax[1].set_xlabel("k"); ax[1].set_ylabel("(-1)^k pi (a_k - A_k)/(sqrt2 tau^2)")
ax[1].set_title("every Fourier coefficient starts at zeta(2) + zeta(3) tau")
ax[1].legend()
plt.tight_layout()
plt.savefig("research/10_frontier/F1_fourier_squircle/f1_boundary_layer.png", dpi=150)
print("figure saved: f1_boundary_layer.png (dotted: uniform formula 1/(4k^2) - pi tau/(2k sinh 2pi k tau), over tau^2)")
