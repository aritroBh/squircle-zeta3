import sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
g=S.init()
def ssub(A,B):
    R=dict(A)
    for k,c in B.items():
        R[k]=R.get(k,Fr(0))-c
        if R[k]==0: del R[k]
    return R
def drop_zero(T):
    return {k:c for k,c in T.items() if all(x!='ZERO' for x in k)}
for (i,j) in [(1,2),(2,5),(6,1),(5,6),(1,6),(6,2)]:
    ri,rj=E.finv(E.X[i]),E.finv(E.X[j])
    s3 = S.S3(E.ZERO, ri, rj)
    s3_nozero = drop_zero(s3)
    sl21,_ = S.sym_Lam(i,j)
    d = ssub(s3_nozero, sl21)
    print(f"Lam({i},{j}): S3(no-zero) terms={len(s3_nozero)}, S_Li21 terms={len(sl21)}, DIFF terms={len(d)}")
    if 0<len(d)<12: print("   ", d)
