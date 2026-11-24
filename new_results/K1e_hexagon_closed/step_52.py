
import mpmath as mp, json, sys, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps = 60
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
from vjet import arc_coeffs_all
import sympy as sp
# verifier polar formulation for the fresh hexagon
V = [sp.Matrix([1,0]), sp.Matrix([3,3]), sp.Matrix([0,1])]
kk = [sp.Rational(1,18)]*3; ww = [sp.Rational(3), sp.Rational(5), sp.Rational(3)]
vs = [(-mp.mpf(int(V[j][1]))/18, mp.mpf(int(V[j][0]))/18) for j in range(3)]
cs = [mp.mpf(int(w)) for w in ww]; lams=[0,0,0]
angs = [mp.atan2(v[1],v[0]) for v in [(1,0),(3,3),(0,1)]] + [mp.pi]
print("angles:", [mp.nstr(a,8) for a in angs])
for q in range(3):
    co = arc_coeffs_all(angs[q], angs[q+1], vs, cs, lams)
    print(f"polar arc {q+1}: P0={mp.nstr(co[0],12)} (L={mp.nstr(mp.sqrt(13) if q<2 else mp.sqrt(2),12)}), P3_arc={mp.nstr(6*co[3],30)}")
