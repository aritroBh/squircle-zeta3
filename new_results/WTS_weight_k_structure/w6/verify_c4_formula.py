"""verify step: formula-vs-stored-data checks at precision HIGHER than the original PSLQ/LLL runs.
 (1) the claimed c_4 closed form (transcribed again from identify_w6.md) vs the stored fitted c_4 in coef_220/250/400/900.json,
     each evaluated at mp.mp.dps = (digits stored + 100);
 (2) the identity  c_4 = (1/3)[y^6]E(y)  between the closed form and the Gamma-function generating function, evaluated at
     mp.mp.dps = 1100 (the original runs used <= 900);
 (3) the generating-function law for weights 2,3 against their exact values.
mp.mp.dps is set AFTER all imports and printed with every residual.
"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp

def E_series(N):
    s = [None] + [(mp.log(2) if n == 1 else (mp.altzeta(n) if n % 2 == 1 else mp.zeta(n))) for n in range(1, N + 1)]
    L = [mp.mpf(0)] + [s[n] / n for n in range(1, N + 1)]
    E = [mp.mpf(1)] + [mp.mpf(0)] * N
    for n in range(1, N + 1):
        E[n] = sum(j * L[j] * E[n - j] for j in range(1, n + 1)) / n
    return E

def c4_closed():
    lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi ** 2
    return (168 * lam ** 6 + 420 * pi2 * lam ** 4 + 378 * pi2 ** 2 * lam ** 2 + 61 * pi2 ** 3 + 5040 * z3 * lam ** 3
            + 2520 * pi2 * z3 * lam + 3780 * z3 ** 2 + 22680 * z5 * lam) / 362880

for run, digits in ((220, 210), (250, 240), (400, 390), (900, 880)):
    d = json.load(open(os.path.join(W, f"coef_{run}.json")))
    mp.mp.dps = digits + 100
    c4s = mp.mpf(d["coefficients"][4])
    cf = c4_closed()
    print(f"run {run}: stored c_4 has {len(d['coefficients'][4].replace('-','').replace('.',''))} sig. digits; "
          f"closed - stored = {mp.nstr(cf - c4s, 5)}   (mp.mp.dps = {mp.mp.dps})", flush=True)

mp.mp.dps = 1100
print("mp.mp.dps =", mp.mp.dps)
E = E_series(8)
cf = c4_closed()
print("closed form - (1/3)[y^6]E(y) =", mp.nstr(cf - E[6] / 3, 5), f"(mp.mp.dps = {mp.mp.dps})")
# weights 2 and 3: g_0 = (1/2)[y^2]E = (pi^2 + 6 L^2)/24 ; g_1 = -(1/2)[y^3]E = -S
lam = mp.log(2); z3 = mp.zeta(3); pi2 = mp.pi ** 2
print("weight 2: (1/2)[y^2]E - (pi^2+6 log^2 2)/24 =", mp.nstr(E[2] / 2 - (pi2 + 6 * lam ** 2) / 24, 5))
S = z3 / 8 + pi2 * lam / 24 + lam ** 3 / 12
print("weight 3: -(1/2)[y^3]E + S              =", mp.nstr(-E[3] / 2 + S, 5))
# Gamma-function form at moderate precision
mp.mp.dps = 60
Ef = lambda y: mp.mpf(2) ** y * mp.gamma(1 - y) * mp.gamma(1 + y / 2) / mp.gamma(1 - y / 2)
tay = mp.taylor(Ef, 0, 8)
E60 = E_series(8)
print(f"Gamma form (mp.taylor) vs series: max |diff| over y^0..y^8 = {mp.nstr(max(abs(tay[n] - E60[n]) for n in range(9)), 5)} (mp.mp.dps = {mp.mp.dps})")
