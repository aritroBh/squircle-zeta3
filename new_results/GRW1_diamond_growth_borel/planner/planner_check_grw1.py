# Planner check of GRW1 (independent of GRW1's code): (1) Lemma 1, I(eta) = int e^{-u} H(eta u) K(eta,u) du with H = sqrt2 sqrt(1+e^{-2z}),
# against the perimeter of |x|^p+|y|^p = 1, p = 1/(1+eta); (2) plain Pade-Borel [m/m] digits from f_0..f_180 (D11b2 data, 450 digits).
import json, mpmath as mp
F = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))["f450"]
def I_perim(eta):
    p = 1 / (1 + eta); s = 2 / p
    f = lambda t: s * mp.sqrt(mp.cos(t)**(2*s - 2) * mp.sin(t)**2 + mp.sin(t)**(2*s - 2) * mp.cos(t)**2)
    P = 8 * mp.quad(f, [0, mp.pi/8, mp.pi/4])
    return P / (4 * mp.sqrt(2) * (1 + eta))
def I_lemma1(eta):
    H = lambda z: mp.sqrt(2) * mp.sqrt(1 + mp.exp(-2*z))
    return mp.quad(lambda u: mp.exp(-u) * H(eta*u) * (1 + mp.exp(-u))**(-2 - eta), [0, 1, 10, 40, mp.inf])
mp.mp.dps = 90
ref = {}
for eta in [mp.mpf('0.1'), mp.mpf('0.2'), mp.mpf('0.3'), mp.mpf('-0.2')]:
    a, b = I_perim(eta), I_lemma1(eta); ref[float(eta)] = a
    print(f"eta={mp.nstr(eta,2)}: I(perimeter) = {mp.nstr(a, 30)}   |perimeter - Lemma 1| = {mp.nstr(abs(a-b), 3)}")
mp.mp.dps = 420
f = [mp.mpf(x) for x in F]; N = len(f) - 1
bcoef = [f[n] / mp.factorial(n) for n in range(N + 1)]
print("f_0..f_%d loaded; f_2 = %s" % (N, mp.nstr(f[2], 20)))
for m in (60, 75, 90):
    P, Q = mp.pade(bcoef[:2*m + 1], m, m)
    mp.mp.dps = 300
    R = lambda z: mp.polyval(P[::-1], z) / mp.polyval(Q[::-1], z)
    row = []
    for eta in (0.1, -0.2):
        e = mp.mpf(str(eta))
        S = mp.quad(lambda t: mp.exp(-t) * R(e * t), [0, 2, 10, 40, 120, 400])
        row.append(f"eta={eta}: {mp.nstr(-mp.log10(abs(S - ref[eta])), 4)} digits")
    print(f"[{m}/{m}] Pade-Borel: " + "; ".join(row))
    mp.mp.dps = 420
