import sys, json, time
sys.path.insert(0, 'research/08_theory/D11b_stokes_numeric')
from mpmath import mp, mpf, mpc, pi, sqrt, quad, exp, nstr
from math import factorial as gfac
from d11b_util import pade_coeffs, pade_eval

t0 = time.time()
mp.dps = 400
f = [mpf(v) for v in json.load(open('research/08_theory/D11b_stokes_numeric/step2_num450_N180.json'))['f450']]
sq2 = sqrt(2)
f1 = [f[n] - sq2*((-1)**n) for n in range(len(f))]
b1 = [f1[n]/mpf(gfac(n)) for n in range(len(f1))]
refs = json.load(open('research/08_theory/D11b_stokes_numeric/refs.json'))

pades = {}
for N in [80, 85, 90]:
    pades[N] = pade_coeffs(b1[:2*N+1], N, N)
print("pade builds: %.1f s" % (time.time()-t0))

def bs1(num, den, eta, dps=100):
    old = mp.dps; mp.dps = dps
    try:
        pts = [0,1,4,16,40,mp.inf] if eta > 0 else list(range(0, 84, 4))+[84]
        v = quad(lambda s: exp(-s)*pade_eval(num, den, eta*s), pts)
    finally:
        mp.dps = old
    return v

print("\n==== PART 1 FINAL: F(eta) = sqrt2/(1+eta) + int_0^inf e^{-s} PB1[90/90](eta s) ds ====")
final = {}
for es in ['0.1', '0.2', '0.3', '-0.2']:
    eta = mpf(es); r = mpf(refs[es])
    v80 = sq2/(1+eta) + bs1(pades[80][0], pades[80][1], eta)
    v85 = sq2/(1+eta) + bs1(pades[85][0], pades[85][1], eta)
    v90 = sq2/(1+eta) + bs1(pades[90][0], pades[90][1], eta)
    res = v90 - r
    digs = -mp.log10(abs(res)/abs(r))
    conv = max(abs(v90-v85), abs(v85-v80))
    final[es] = (v90, res, digs, conv)
    print(f"\neta = {es}")
    print(f"  F_PadeBorel = {nstr(v90, 60)}")
    print(f"  F_direct    = {nstr(r, 60)}")
    print(f"  residual = {nstr(res,3)}  rel = {nstr(abs(res)/abs(r),3)}  digits = {nstr(digs,5)}")
    print(f"  order-convergence |v90-v85|,|v85-v80| max = {nstr(conv,3)}")
num, den = pades[90]
json.dump({'num': [nstr(v, 350) for v in num], 'den': [nstr(v, 350) for v in den]},
          open('research/08_theory/D11b_stokes_numeric/pade1_9090.json','w'))
json.dump({k: nstr(v[0], 70) for k, v in final.items()},
          open('research/08_theory/D11b_stokes_numeric/part1_final.json','w'))
print("\nsaved pade1_9090.json, part1_final.json;  total %.1f s" % (time.time()-t0))
