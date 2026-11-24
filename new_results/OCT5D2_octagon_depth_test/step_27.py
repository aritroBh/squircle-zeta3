import os, sys
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from qsqrt2 import QS2 as Q
Z = Q(0); ONE_ = Q(1); P = Q(2,1); Qq = Q(-1,-1)

def integrability_residuals(s):
    """For each adjacent position i and each (n-2)-context u, sum s(w)*det(w_i,w_{i+1}).
    Must be 0 for an integrable (iterated-integral) symbol."""
    from collections import defaultdict
    bad = []
    n = len(next(iter(s)))
    for i in range(n-1):
        groups = defaultdict(lambda: Z)
        for w, c in s.items():
            ctx = w[:i]+w[i+2:]
            det = (w[i+1]-w[i])  # det(e_a,e_b) = b - a for basis X=0,Y=1 in 2D: det = a*1... 
            # det(v_i, v_{i+1}) with v in {X=(1,0), Y=(0,1)}: det = a_i*b_{i+1}-a_{i+1}*b_i = w[i+1]-w[i]
            groups[ctx] += c*det
        for ctx, v in groups.items():
            if v != Z: bad.append((i, ctx, v))
    return bad

for wtest in [(P,Qq,Z,ONE_,P), (Qq,P,P,Z,ONE_), (Z,P,Qq,ONE_,ONE_)]:
    s5 = S.gsym(wtest)
    bad = integrability_residuals(s5)
    print(wtest, "words:", len(s5), "integrability violations:", len(bad))
    for b in bad[:5]: print("   ", b)
