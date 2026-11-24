
import sys; sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr

S22 = S_Li22(); S31 = S_Li31(); S13 = S_Li13()
Dm = delta22(specialize_lin(S22, spec_map(-1)))   # y=1-c
Dp = delta22(specialize_lin(S22, spec_map(+1)))   # y=1+c
D31m = delta22(specialize_lin(S31, spec_map(-1)))
D31p = delta22(specialize_lin(S31, spec_map(+1)))
print("d22 Li31(1-c) == -d22 Li22(1-c)?", D31m == d22_scale(Dm,-1))
print("d22 Li31(1+c) == -d22 Li22(1+c)?", D31p == d22_scale(Dp,-1))

# DT_A = ((1-c)/c) Li22(1/2,1-c) - ((1+c)/c) Li22(1/2,1+c)
# coefficients as Laurent polys in c: (1-c)/c = 1/c - 1 ; (1+c)/c = 1/c + 1
def lau_scale(d, poly):  # poly: dict e->Fraction
    out = {}
    for k,v in d.items():
        out[k] = {e: v*c0 for e,c0 in poly.items()}
    return out
def lau_add(d1, d2):
    out = dict(d1)
    for k,v in d2.items():
        cur = out.setdefault(k, {})
        for e,c0 in v.items():
            cur[e] = cur.get(e, Fr(0)) + c0
            if cur[e]==0: del cur[e]
        if not cur: del out[k]
    return out
A_minus = lau_scale(Dm, {-1: Fr(1), 0: Fr(-1)})     # (1-c)/c
A_plus  = lau_scale(Dp, {-1: Fr(-1), 0: Fr(-1)})   # -(1+c)/c
DTA = lau_add(A_minus, A_plus)
# DT_B = ((1-c)/c) Li22(1/2,1-c) + ((1+c)/c) Li31(1/2,1+c)
B_plus = lau_scale(D31p, {-1: Fr(1), 0: Fr(1)})
DTB = lau_add(A_minus, B_plus)
print("DT_A == DT_B (with [2] kept)?", DTA == DTB)

# target
T = lau_add(lau_scale({(3,5): Fr(1)}, {0: Fr(1), -1: Fr(-1)}), lau_scale({(4,5): Fr(1)}, {0: Fr(-1), -1: Fr(-1)}))
resid = lau_add(DTA, lau_scale(T, {0: Fr(-1)}))
print("\nResidual d22(DT_A) - d22(Omega), with [2] kept:")
for k,v in sorted(resid.items()):
    print(f"  {WNAME[k[0]]}^{WNAME[k[1]]}: ", {e: str(c) for e,c in v.items()})
resid0 = {k: v for k,v in resid.items() if any(i in (3,4,5) for i in k)}  # components not involving m=[2]
# check all residual components involve index 0 (letter m)
allm = all(0 in k for k in resid)
print("all residual components involve [2]?", allm)
print("\nwith [2] dropped: residual components:",
      {k:v for k,v in resid.items() if 0 not in k})
open('progress.log','a').write("step_05 DT_A=DT_B exact; residual vs Omega = 4 [2]-components; [2]-dropped exact +sign\n")
