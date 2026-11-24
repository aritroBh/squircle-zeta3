
import pickle, sys
sys.path.insert(0,'.')
from fractions import Fraction as Fr
s2 = pickle.load(open('../COR8_depth_two_any_angle/symbols_2kept.pkl','rb'))
DT_cor8 = s2['DT']
mine = pickle.load(open('syms_omega_dt.pkl','rb'))
S_DTA = mine['S_DTA']
# normalize keys: mine are tuples of tuples of ints; cor8 tuples of Fractions
def norm(sym):
    out = {}
    for k,v in sym.items():
        nk = tuple(tuple(int(x) for x in vv) for vv in k)
        out[nk] = {e: Fr(c0) for e,c0 in v.items()}
    return out
A = norm(DT_cor8); B = norm(S_DTA)
keysA, keysB = set(A), set(B)
print("cor8 DT tensors:", len(A), " mine:", len(B))
print("common:", len(keysA & keysB))
print("\nOnly in COR8:")
for k in sorted(keysA - keysB):
    print("   ", k, {e:str(c0) for e,c0 in A[k].items()})
print("\nOnly in mine:")
for k in sorted(keysB - keysA):
    print("   ", k, {e:str(c0) for e,c0 in B[k].items()})
print("\nDifferent coefficients:")
for k in sorted(keysA & keysB):
    if A[k] != B[k]:
        print("   ", k, " cor8:", {e:str(c0) for e,c0 in A[k].items()}, " mine:", {e:str(c0) for e,c0 in B[k].items()})
