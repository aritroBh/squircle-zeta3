
import sys, time
sys.path.insert(0, '<project>/research/08_theory/D16_pure_weight')
import sympy as sp
from d16_engine import full_S, s_log, s_mul, Lambda, logRat, build_B

t0 = time.time()
n, K = 3, 12
B, eLam, S = full_S(n, K)          # S = S_n/S_n(1)/Rat_n in tau, formal Z's
W = s_log(S, K)                    # w_k: log[(2-p)S/(4sqrt3)] = sum w_k tau^k
print(f"exact engine n=3 K=12: {time.time()-t0:.1f}s")
# numeric evaluation
import mpmath as mp
mp.mp.dps = 40
def ev(d):
    out = {}
    for k, md in d.items():
        v = mp.mpf(0)
        for m, c in md.items():
            t = mp.mpf(sp.Rational(c))
            for i in m: t *= mp.zeta(i)
            v += t
        out[k] = v
    return out
Wn = ev(W)
paper = {2: -mp.pi**2/9, 3: -16*mp.zeta(3)/9, 4: -37*mp.pi**4/810, 5: -160*mp.zeta(5)/27,
         6: 16*mp.zeta(3)**2/81 - 572*mp.pi**6/76545}
print("w_k (log series) vs paper thm:oct:")
for k in range(2, 7):
    print(f"  k={k}: engine={mp.nstr(Wn[k],25)}  paper={mp.nstr(paper[k],25)}  diff={mp.nstr(abs(Wn[k]-paper[k]),3)}")
for k in range(7, 13):
    print(f"  k={k}: {mp.nstr(Wn[k],25)}")
# also save exact W
import json
json.dump({str(k): {' '.join(map(str,m)): str(c) for m,c in W[k].items()} for k in W},
          open('w_exact_n3_K12.json','w'), indent=1)
print("saved w_exact_n3_K12.json")
