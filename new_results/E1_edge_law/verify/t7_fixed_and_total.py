# T7: (a) the same rational hexagon (b=0) in the TASK's fixed-weight convention (unit normals, lam = 0):
#        K_p = { (1/4)|y|^p + (3 sqrt2/4)|(y-x)/sqrt2|^p + (1/4)|x|^p <= 1 }   (same K_1, different p-family)
#     per-arc PSLQ at two precisions;
#     (b) PSLQ on the TOTAL P'''(1) of the hexagon (power-type family of t3) at two precisions, basis {sqrt5,sqrt2} x 7.
import sys, json, mpmath as mp
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
from vjet import *
out = {}
for dps in (70, 110):
    mp.mp.dps = dps
    z3, l2, pi2 = mp.zeta(3), mp.log(2), mp.pi**2
    s2, s5 = mp.sqrt(2), mp.sqrt(5)
    vs = [(mp.mpf(0), mp.mpf(1)), (-1/s2, 1/s2), (mp.mpf(1), mp.mpf(0))]
    cs = [mp.mpf(1)/4, 3*s2/4, mp.mpf(1)/4]; lams = [0, 0, 0]
    angs = [mp.mpf(0), mp.pi/4, mp.pi/2, mp.pi]
    Ls = [s5, s5, s2]; Z = [18*s5/125, 18*s5/125, -15*s2/128]
    print(f"== (a) fixed-weight convention, dps {dps}")
    tot = mp.mpf(0)
    for q in range(3):
        co = arc_coeffs_all(angs[q], angs[q+1], vs, cs, lams)
        P3 = 6*co[3]; tot += 2*P3
        rel = mp.pslq([P3/Ls[q], z3, pi2*l2, l2**3, pi2, l2**2, l2, 1], maxcoeff=10**12, maxsteps=10**6)
        cz = -mp.mpf(rel[1])/rel[0]*Ls[q] if rel else None
        print(f"   arc {q+1}: P'''_arc = {mp.nstr(P3, 30)}  PSLQ {rel}  z3-coef = {mp.nstr(cz, 20) if rel else None}  edge law {mp.nstr(Z[q], 20)}")
        out[f'fixed_{dps}_arc{q+1}'] = dict(P3=mp.nstr(P3, dps), rel=rel)
    print(f"   total P'''(1) (fixed-weight family) = {mp.nstr(tot, 40)}")
    # (b) total of the t3 family from its saved json
    d = json.load(open('research/10_frontier/E1_edge_law/verify/out/t3_b0_110.json'))
    Pt = mp.mpf(d['P3_total'])
    B = [z3, pi2*l2, l2**3, pi2, l2**2, l2, mp.mpf(1)]
    rel = mp.pslq([Pt] + [s5*b for b in B] + [s2*b for b in B], maxcoeff=10**12, maxsteps=10**7)
    print(f"== (b) TOTAL P'''(1) of the t3 hexagon (b=0) at {dps} digits: PSLQ on [P''', sqrt5*(z3,pi2 log2,log2^3,pi2,log2^2,log2,1), sqrt2*(...)]: {rel}")
    if rel:
        print(f"   total zeta(3) coefficient = ({-rel[1]}/{rel[0]}) sqrt5 + ({-rel[8]}/{rel[0]}) sqrt2 ; edge law sum Z_e = (72/125) sqrt5 - (15/64) sqrt2")
    Rt = Pt - (72*s5/125 - 15*s2/64)*z3
    rel2 = mp.pslq([Rt] + [s5*b for b in B[1:]] + [s2*b for b in B[1:]], maxcoeff=10**12, maxsteps=10**7)
    print(f"   P''' - (sum Z_e) z3 on the z3-free basis {{sqrt5,sqrt2}} x 6: {rel2}")
    out[f'total_{dps}'] = dict(rel=rel, rel_without_z3=rel2)
json.dump(out, open('research/10_frontier/E1_edge_law/verify/out/t7.json', 'w'), indent=1, default=str)
