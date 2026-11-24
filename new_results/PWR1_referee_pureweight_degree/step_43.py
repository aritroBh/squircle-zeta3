import mpmath as mp, json, sympy as sp
mp.mp.dps = 100
import pwr1_engine as E

w = json.load(open('w_exact_n2_tau10.json'))
znum = {E.Z[k]: mp.zeta(k) for k in range(2, 12)}
W = {int(k): mp.mpf(sp.N(sp.sympify(v), 50).subs({s: m for s, m in znum.items()}) if False else sp.sympify(v).subs(znum).evalf(60)) for k, v in w.items()}
vals = {mp.mpf(k): mp.mpf(v) for k, v in json.load(open('quad_logP_n2_hi.json')).items()}

def partial(tau, kmax):
    return mp.fsum(W[k]*tau**k for k in range(2, kmax+1))

print("remainder test: R(tau) = logP - sum_{2..8} w_k tau^k; ratio R/(w9 tau^9) -> 1 ?")
for tau, v in sorted(vals.items()):
    R = v - partial(tau, 8)
    pred = W[9]*tau**9
    print(f"tau={mp.nstr(tau,3)}: R={mp.nstr(R,12)}  R/(w9 tau^9)={mp.nstr(R/pred, 12)}")

print("\nnext-order check: (R - w9 tau^9)/(w10 tau^10) -> 1 ?")
for tau, v in sorted(vals.items()):
    R = v - partial(tau, 9)
    pred = W[10]*tau**10
    print(f"tau={mp.nstr(tau,3)}: ratio={mp.nstr(R/pred, 12)}")
