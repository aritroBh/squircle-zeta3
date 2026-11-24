import sys
sys.path.insert(0,"<project>/research/08_theory/CHN1_chain_bodies")
from fractions import Fraction as F
import engine
def smap(a):
    a1,a2,a3,a4=a; return (a2,F(0),a3,a2+a4,a1+a3)
def eng(a): return engine.reduced_coeffs(tuple(F(x) for x in a),7)
def M_12(a,W):
    a=tuple(F(x) for x in a); s=smap(a); e=eng(a)
    # single-zeta coefficient index 0 of BASIS[W] is zeta_W (odd W)
    m11=e[W][0]
    shift=((s[3]+s[4])**W-s[3]**W-s[4]**W)/W    # 4pt Gamma-ratio coefficient of zeta_W at s1'=0 (odd W)
    return m11-shift
for W in (3,5,7):
    for a in [(1,1,2,2),(2,2,1,1),(1,1,3,3),(3,3,1,1),(1,2,1,2),(2,1,2,1),(1,2,2,1)]:
        perm=(a[2],a[3],a[0],a[1])
        print("W=%d a=%s: M_12(a)=%s  M_21(a)=M_12(perm a)=%s  equal: %s"%(W,a,M_12(a,W),M_12(perm,W),M_12(a,W)==M_12(perm,W)))
