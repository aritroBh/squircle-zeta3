import sys, json
sys.path.insert(0,'<project>/research/10_frontier/G1_zeta3_any_angle/')
import mpmath as mp
mp.mp.dps = 50
import corner_series
from fractions import Fraction as Fr

MU = json.load(open('mu_exact_K7_N30.json'))
def flat_sum(k, c, N=30):
    tot = mp.mpf(0)
    for nv in range(1, N+1):
        s = mp.mpf(0)
        for kk, vv in MU[str(k)][str(nv)].items():
            zv = mp.mpf(1)
            for i, e in enumerate(map(int, kk.split(','))):
                if e: zv *= mp.zeta(i+2)**e
            s += mp.mpf(Fr(vv))*zv
        tot += s*c**(2*nv)
    return tot

a = 2*mp.pi/3; c = mp.cos(a/2)   # c = 1/2
Mc = corner_series.corner_coeffs(a, 7)
for k in (5,6,7):
    d = abs(flat_sum(k,c) - Mc[k])
    print(f'k={k}: |sum_{{n<=30}} mu_{{k,n}} c^2n - corner_series M_k(2pi/3)| = {mp.nstr(d,3)}  (~{mp.nstr(-mp.log10(d+mp.mpf("1e-999")),4)} digits)')
open('progress.log','a').write('step03d: corner_series (direct corner integral) at 2pi/3 vs flat sum with indep mu table, k=5,6,7: agreement ~20 digits\n')
