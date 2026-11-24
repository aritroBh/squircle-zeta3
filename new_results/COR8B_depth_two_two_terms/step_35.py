
import sys; sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr
S22 = S_Li22()
Dm = delta22(specialize_lin(S22, spec_map(-1)))
Dp = delta22(specialize_lin(S22, spec_map(+1)))
def lau_scale(d, poly):
    out = {}
    for k,v in d.items():
        if isinstance(v, dict):
            cur = {}
            for e0,cv in v.items():
                for e,c0 in poly.items():
                    cur[e0+e] = cur.get(e0+e, Fr(0)) + cv*c0
            out[k] = {e:c for e,c in cur.items() if c}
        else:
            out[k] = {e: v*c0 for e,c0 in poly.items() if v*c0}
    return {k:v for k,v in out.items() if v}
def lau_add(d1, d2):
    out = {k: dict(v) for k,v in d1.items()}
    for k,v in d2.items():
        cur = out.setdefault(k, {})
        for e,c0 in v.items():
            cur[e] = cur.get(e, Fr(0)) + c0
            if cur[e]==0: del cur[e]
        if not cur: del out[k]
    return out
A_minus = lau_scale(Dm, {-1: Fr(1), 0: Fr(-1)})
A_plus  = lau_scale(Dp, {-1: Fr(-1), 0: Fr(-1)})
DTA = lau_add(A_minus, A_plus)
T = {(3,5): {0: Fr(1), -1: Fr(-1)}, (4,5): {0: Fr(-1), -1: Fr(-1)}}
resid = lau_add(DTA, lau_scale(T, {0: Fr(-1)}))
print("Residual d22(DT_A) - d22(Omega), with [2] kept:")
for k,v in sorted(resid.items()):
    print(f"  ({WNAME[k[0]]})^({WNAME[k[1]]}): ", {e: str(c) for e,c in v.items()})
print("all residual components involve letter [2] (index 0)?", all(0 in k for k in resid))
print("number of residual components:", len(resid))
json = __import__('json')
json.dump({f"{k[0]},{k[1]}": {str(e): str(c) for e,c in v.items()} for k,v in resid.items()},
          open('DT_A_residual_2kept.json','w'), indent=1)
print("saved DT_A_residual_2kept.json")
