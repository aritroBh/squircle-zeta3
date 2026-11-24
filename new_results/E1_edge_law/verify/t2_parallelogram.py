# T2: verifier's fully-local test polygon (non-rhombus parallelogram, FIXED weights, not the power convention):
#   K_p = { (|x|^p + |x+2y|^p)/2 <= 1 },  vertices (2,-1),(0,1),(-2,1),(0,-1).
# Per-arc P'''(1) by the verifier's jet code; PSLQ of P'''_e/L_e on {z3, pi^2 log2, log^3 2, pi^2, log^2 2, log2, 1}.
import sys, json, mpmath as mp
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
from vjet import *
dps = int(sys.argv[1]); mp.mp.dps = dps
vs = [(mp.mpf(1), mp.mpf(0)), (mp.mpf(1), mp.mpf(2))]; cs = [mp.mpf(1)/2, mp.mpf(1)/2]; lams = [0, 0]
z3, l2, pi2 = mp.zeta(3), mp.log(2), mp.pi**2
angs = [mp.atan2(-1, 2), mp.pi/2, mp.atan2(1, -2)]
Ls = [2*mp.sqrt(2), mp.mpf(2)]
Zpred = [6*mp.mpf(1)/2/(2*mp.sqrt(2)), 6*mp.mpf(1)/2]   # 6 h^2/L : (h,L) = (1/sqrt2, 2sqrt2), (1, 2)
out = dict(dps=dps, arcs=[])
tot = mp.mpf(0)
for q in range(2):
    co = arc_coeffs_all(angs[q], angs[q+1], vs, cs, lams)
    P3 = 6*co[3]; tot += 2*P3
    X = P3/Ls[q]
    rel = mp.pslq([X, z3, pi2*l2, l2**3, pi2, l2**2, l2, 1], maxcoeff=10**12, maxsteps=10**6)
    cz3 = -mp.mpf(rel[1])/rel[0]*Ls[q] if rel else None
    print(f"arc {q+1}: P0 = {mp.nstr(co[0], 20)} (L = {mp.nstr(Ls[q], 20)}), P'''_arc = {mp.nstr(P3, 35)}")
    print(f"   PSLQ [X, z3, pi2 log2, log2^3, pi2, log2^2, log2, 1]: {rel}")
    if rel:
        print(f"   zeta(3) coefficient = {mp.nstr(cz3, 25)} ; edge law 6h^2/L = {mp.nstr(Zpred[q], 25)} ; diff = {mp.nstr(cz3 - Zpred[q], 3)}")
    out['arcs'].append(dict(P3=mp.nstr(P3, dps), rel=rel, z3=mp.nstr(cz3, 30) if rel else None, Z=mp.nstr(Zpred[q], 30)))
print("total P'''(1) =", mp.nstr(tot, 40), "; edge-law total zeta(3) coeff =", mp.nstr(2*sum(Zpred), 25), "= 3sqrt2/2+6 :", mp.nstr(3*mp.sqrt(2)/2 + 6, 25))
out['P3_total'] = mp.nstr(tot, dps)
json.dump(out, open(f'research/10_frontier/E1_edge_law/verify/out/t2_{dps}.json', 'w'), indent=1, default=str)
