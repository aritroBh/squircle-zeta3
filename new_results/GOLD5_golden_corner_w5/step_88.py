import sys, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
rho  = (Fr(-1,2), Fr(1,2))    # 1/phi
mrho = (Fr(1,2), Fr(-1,2))    # -1/phi
rho3 = (Fr(-2), Fr(1))        # sqrt5-2
rho6 = (Fr(9), Fr(-4))        # 9-4sqrt5
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfscale(q, a):
    return (q*a[0], q*a[1])
# ZQ5 Li5 part of T = L(chi5,5)/sqrt5:
combo = [(Fr(1728, 6875), rho), (Fr(-1728, 6875), mrho), (Fr(-512, 6875), rho3), (Fr(16, 6875), rho6)]
tot = [(Fr(0), Fr(0))]*6
for c, z in combo:
    cl = li5cl.get(z)
    print(f"class Li5({z}): {[G.nf_str(q) for q in cl] if cl else 'MISSING'}")
    if cl:
        tot = [tuple(a + b for a, b in zip(t, nfscale(c, q))) for t, q in zip(tot, cl)]
print("ZQ5 Li5-combo class:", [G.nf_str(q) for q in tot])
print("  zero?", all(not any(q) for q in tot))
