# s07: independent isolation of zeta(3) for the non-regular hexagon H_p = {|x|^p+|y|^p+|x-y|^p <= 1}
# (and its fixed-weight version) by PSLQ on P'''(1) computed with D15b's POLAR coefficient integrals
# (poly_coeffs.P_coeffs_general; no use of the radial-projection derivation), at two precisions.
import sys, json, mpmath as mp
sys.path.insert(0, 'research/08_theory/D15b_polygon_law')
import poly_coeffs as pc
res = {}
for dps in (int(sys.argv[1]), int(sys.argv[2])):
    mp.mp.dps = dps + 20
    r2 = mp.sqrt(2); l2 = mp.log(2)
    for name, lams in (('H_power', [0, 0, mp.log(r2)]), ('H_fixed', [0, 0, 0])):
        bs = [mp.mpf(0), mp.pi/2, mp.pi/4]; ws = [1, 1, r2]
        P3 = 6*pc.P_coeffs_general(bs, ws, lams, kmax=3)[3]
        mp.mp.dps = dps
        z3 = mp.zeta(3); pi2 = mp.pi**2
        names = ['P3', 'z3', 'r2*z3', 'pi2', 'r2*pi2', 'pi2*l2', 'r2*pi2*l2', 'l2^3', 'r2*l2^3',
                 'l2^2', 'r2*l2^2', 'l2', 'r2*l2', '1', 'r2']
        vec = [P3, z3, r2*z3, pi2, r2*pi2, pi2*l2, r2*pi2*l2, l2**3, r2*l2**3, l2**2, r2*l2**2, l2, r2*l2, 1, r2]
        rel = mp.pslq(vec, maxcoeff=10**6, maxsteps=10**6)
        resid = mp.fsum([a*b for a, b in zip(rel, vec)]) if rel else None
        print(f"dps {dps} {name}: P'''(1) = {mp.nstr(P3, 40)}")
        print("   relation:", dict(zip(names, rel)) if rel else None, " residual", mp.nstr(resid, 3) if rel else '')
        if rel:
            c = -mp.mpf(rel[1])/rel[0], -mp.mpf(rel[2])/rel[0]
            print(f"   => zeta(3) coefficient = {mp.nstr(c[0], 15)} + {mp.nstr(c[1], 15)}*sqrt2 ;"
                  f" edge law -3/2 - 3sqrt2/16 = {mp.nstr(-mp.mpf(3)/2, 15)} + {mp.nstr(-mp.mpf(3)/16, 15)}*sqrt2")
        res[f'{name}_{dps}'] = dict(P3=str(P3), rel=rel)
        mp.mp.dps = dps + 20
json.dump(res, open('s07_pslq_H.json', 'w'), indent=1)
