# T1: sanity of the verifier's jet code against the solver's exact closed forms (regular hexagon, H_p power, rhombus)
import sys, time, mpmath as mp
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
from vjet import *
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mp.mp.dps = dps
z2, z3, l2 = mp.zeta(2), mp.zeta(3), mp.log(2)
s3, s2 = mp.sqrt(3), mp.sqrt(2)
def total(vs, cs, lams):
    ks = kink_angles(vs)
    pts = ks + [ks[0] + 2*mp.pi]
    tot = [mp.mpf(0)]*4
    for a, b in zip(pts[:-1], pts[1:]):
        c = arc_coeffs_all(a, b, vs, cs, lams)
        tot = [x+y for x, y in zip(tot, c)]
    return [tot[0], tot[1], 2*tot[2], 6*tot[3]]
cases = []
# regular hexagon, D15: n_i = e(i pi/3)
vs = [(mp.cos(i*mp.pi/3), mp.sin(i*mp.pi/3)) for i in (1, 2, 3)]
exact = z2*(-9*s3/2 + 9*s3*l2/4) - 9*s3/8*z3 - 12*s3*l2**2 + 2*s3*l2**3 + 12*s3*l2
cases.append(("regular hexagon", vs, [1,1,1], [0,0,0], [2*s3, 2*s3*l2, 3*s3*z2/4 - 4*s3*l2 + 2*s3*l2**2, exact]))
# H_p power |x|^p+|y|^p+|x-y|^p
vs = [(1, 0), (0, 1), (1, -1)]
exact = z2*(-6 - 3*s2/4 + 3*s2*l2/8 + 3*l2) + z3*(-mp.mpf(3)/2 - 3*s2/16) - 12*l2**2 - 6*s2*l2**2 + s2*l2**3 + 2*l2**3 + 6*s2*l2 + 12*l2
cases.append(("H_p power", vs, [1,1,1], [0,0,0], [2+s2, None, None, exact]))
# rhombus |x|^p + |y/2|^p (R2: P'''=24a^2b^2/r^3 (z3 - pi^2/3), a=1,b=2)
vs = [(1, 0), (0, mp.mpf(1)/2)]
r = mp.sqrt(5)
cases.append(("rhombus a=1 b=2", vs, [1,1], [0,0], [4*r, 0, None, 24*4/r**3*(z3 - mp.pi**2/3)]))
for name, vs, cs, lams, ex in cases:
    t0 = time.time()
    P = total(vs, cs, lams)
    print(f"== {name} (dps {dps}, {time.time()-t0:.1f}s)")
    for k in range(4):
        if ex[k] is not None:
            print(f"   P^({k})(1) = {mp.nstr(P[k], 30)}   |diff vs solver closed form| = {mp.nstr(abs(P[k]-ex[k]), 3)}")
        else:
            print(f"   P^({k})(1) = {mp.nstr(P[k], 30)}")
