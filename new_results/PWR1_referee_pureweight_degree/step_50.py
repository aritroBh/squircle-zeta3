import mpmath as mp, json, sympy as sp
mp.mp.dps = 60
import pwr1_engine as E
w3 = json.load(open('w_exact_n3_tau10.json'))
znum = {E.Z[k]: mp.zeta(k) for k in range(2, 12)}
W3 = {int(k): mp.mpf(str(sp.sympify(v).subs(znum).evalf(50))) for k, v in w3.items()}
vals3 = {mp.mpf(k): mp.mpf(v) for k, v in json.load(open('quad_logS_n3.json')).items()}

def logRat3(tau):
    return mp.log(1-tau) - mp.log(1-2*tau)   # Rat_3 = (1-tau)/(1-2 tau) = 1/(2-p)

print("n=3: pure part from quadrature = logS - logRat3; remainder after sum_{2..8}:")
for tau, v in sorted(vals3.items()):
    pure_num = v - logRat3(tau)
    R = pure_num - mp.fsum(W3[k]*tau**k for k in range(2, 9))
    pred = W3[9]*tau**9
    print(f"tau={mp.nstr(tau,3)}: R={mp.nstr(R,10)}  R/(w9 tau^9)={mp.nstr(R/pred, 10)}")

print("\ncoefficient extraction w2..w5 (successive stripping, tau=0.001):")
tau = mp.mpf('0.001'); v = vals3[tau] - logRat3(tau)
rem = v
for k in range(2, 6):
    wk = rem/tau**k
    err = abs(wk - W3[k]); rel = err/abs(W3[k])
    print(f"  w_{k}(3): extracted={mp.nstr(wk,25)} exact={mp.nstr(W3[k],25)} relerr={mp.nstr(rel,2)} ({mp.nstr(-mp.log10(rel),3)} digits)")
    rem = rem - W3[k]*tau**k
open('progress.log','a').write("step_46-50: n=2 quadrature 100dps at 9 taus + 70dps at 7 taus; remainder test R/(w9 tau^9)->1, (R-w9t9)/(w10t10)->1 PASS; LSQ fit w2..w12: 20.5..0.5 digits; purity residuals w2 1e-21, w4 8e-15, w5 7e-12; n=3 engine tau^10 exact, gamma-free pure weight, matches H1 T5 L2..L7 exactly\n")
