# T3: verifier's rational hexagon, NOT fully local, per-arc P'''(1) by the verifier's own jet code,
# then PSLQ of P'''_e / L_e on the weight<=3 alternating-MZV basis {z3, pi^2 log2, log^3 2, pi^2, log^2 2, log2, 1}.
# Hexagon (shear parameter b): V1=(1,0), V2=(2+2b,2), V3=(b,1).  N(x) = sum_j c_j |cross(V_j,x)|, c=(1/4,3/8,1/4).
# p-family: F_p = |cross(V1,x)/4|^p + 3*|cross(V2,x)/8 *... | -- written below as c_j' |<v_j,x>|^p with rational v_j
# chosen so that every |<v_j,V_k>| is 1/4 or 1/2 (only log 2 enters the mu's).
import sys, json, time, mpmath as mp, sympy as sp
sys.path.insert(0, 'research/10_frontier/E1_edge_law/verify')
from vjet import *
b = sp.Rational(sys.argv[1]); dps = int(sys.argv[2])
mp.mp.dps = dps
V = [sp.Matrix([1, 0]), sp.Matrix([2+2*b, 2]), sp.Matrix([b, 1])]
cr = lambda P, Q: P[0]*Q[1] - P[1]*Q[0]
# <v_j, x> := cross(V_j, x) * k_j ; F_1 = sum_j cw_j |<v_j,x>|
c = [sp.Rational(1, 4), sp.Rational(3, 8), sp.Rational(1, 4)]
# middle term: (3/8)|cross(V2,x)| = 3 * |cross(V2,x)/8|  -> v_2 = J V2 / 8, weight 3 ; others weight 1, v_j = J V_j / 4
kk = [sp.Rational(1, 4), sp.Rational(1, 8), sp.Rational(1, 4)]; ww = [1, 3, 1]
vs_ex = [(-V[j][1]*kk[j], V[j][0]*kk[j]) for j in range(3)]     # <v,x> = kk*cross(V_j,x)
for j in range(3):
    for q in range(3):
        val = abs(kk[j]*cr(V[j], V[q]))
        assert val in (0, sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4), 1), val
# check N(V_q) = 1
Nf = lambda x: sum(ww[j]*abs(vs_ex[j][0]*x[0] + vs_ex[j][1]*x[1]) for j in range(3))
assert all(Nf(Vq) == 1 for Vq in V), [Nf(Vq) for Vq in V]
print("vertices:", [tuple(v) for v in V], " |<v_j,V_k>| values:", sorted(set(abs(kk[j]*cr(V[j], V[q])) for j in range(3) for q in range(3))))
# exact edge-law data (verifier's own implementation of the boxed formula AND of the geometric form)
verts = V + [-v for v in V]
def edge_law(V0, V1):
    D = V1 - V0; L = sp.sqrt(D.dot(D)); tau = D/L
    h = abs(cr(V0, tau)); s0 = V0.dot(tau); u = -s0/L; v = 1 - u
    j0 = [j for j in range(3) if cr(V[j], V0) == 0][0]; j1 = [j for j in range(3) if cr(V[j], V1) == 0][0]
    al = ww[j0]*abs(kk[j0]*cr(V[j0], V1)); be = ww[j1]*abs(kk[j1]*cr(V[j1], V0))
    Z = sp.nsimplify(sp.simplify(6*al*be*h**2/L*(3*(al*u + be*v) - 2)))
    other = [j for j in range(3) if j not in (j0, j1)][0]
    r = sp.simplify(cr(V[other], V1)/cr(V[other], V0))
    return dict(L=L, h=h, u=u, al=al, be=be, Z=Z, r=r)
eds = [edge_law(verts[q], verts[(q+1) % 6]) for q in range(3)]
for q, e in enumerate(eds):
    print(f"edge {q+1}: L={e['L']}, h={sp.simplify(e['h'])}, u={sp.simplify(e['u'])}, alpha={e['al']}, beta={e['be']}, r={e['r']}, Z_e={e['Z']} = {sp.N(e['Z'], 25)}")
Ztot = sp.nsimplify(sp.simplify(2*sum(e['Z'] for e in eds)))
print("edge-law total sum_e Z_e (6 edges) =", Ztot, "=", sp.N(Ztot, 30))
# numerics: per arc
vs = [(mp.mpf(sp.Rational(x).p)/sp.Rational(x).q, mp.mpf(sp.Rational(y).p)/sp.Rational(y).q) for x, y in vs_ex]
cs = [mp.mpf(w) for w in ww]; lams = [0, 0, 0]
angs = [mp.atan2(mp.mpf(sp.Rational(v[1]).p)/sp.Rational(v[1]).q, mp.mpf(sp.Rational(v[0]).p)/sp.Rational(v[0]).q) for v in V]
angs.append(mp.pi)
z3, l2, pi2 = mp.zeta(3), mp.log(2), mp.pi**2
basis_names = ['X', 'z3', 'pi2*log2', 'log2^3', 'pi2', 'log2^2', 'log2', '1']
out = dict(b=str(b), dps=dps, arcs=[])
tot3 = mp.mpf(0)
for q in range(3):
    t0 = time.time()
    co = arc_coeffs_all(angs[q], angs[q+1], vs, cs, lams)
    P3 = 6*co[3]
    tot3 += 2*P3
    Le = mp.sqrt(mp.mpf(sp.Rational(eds[q]['L']**2)))
    X = P3/Le
    rel = mp.pslq([X, z3, pi2*l2, l2**3, pi2, l2**2, l2, 1], maxcoeff=10**15, maxsteps=10**6)
    Zn = mp.mpf(sp.N(eds[q]['Z'], dps+10))
    info = dict(edge=q+1, P0=mp.nstr(co[0], dps), P3_arc=mp.nstr(P3, dps), relation=rel, secs=time.time()-t0)
    print(f"-- arc {q+1} ({time.time()-t0:.1f}s): P0_arc = {mp.nstr(co[0], 25)} (exact L_e = {mp.nstr(Le, 25)}), P'''_arc = {mp.nstr(P3, 40)}")
    if rel:
        cz3 = -mp.mpf(rel[1])/rel[0]*Le
        print(f"   PSLQ relation on {basis_names}: {rel}")
        print(f"   => zeta(3) coefficient of P'''_arc = {sp.nsimplify(sp.Rational(-rel[1], rel[0]))} * L_e = {mp.nstr(cz3, 25)};  edge-law Z_e = {mp.nstr(Zn, 25)};  difference = {mp.nstr(cz3 - Zn, 25)}")
        resid = mp.fsum(r*x for r, x in zip(rel, [X, z3, pi2*l2, l2**3, pi2, l2**2, l2, 1]))
        print(f"   residual of relation: {mp.nstr(resid, 3)}")
        info['z3coef_pslq'] = mp.nstr(cz3, 30); info['Z_e'] = mp.nstr(Zn, 30); info['diff'] = mp.nstr(cz3 - Zn, 30)
    else:
        print("   PSLQ: no relation found")
    # zeta3-free test of the remainder P'''_arc - Z_e*zeta3 on the basis WITHOUT zeta(3)
    R = (P3 - Zn*z3)/Le
    rel2 = mp.pslq([R, pi2*l2, l2**3, pi2, l2**2, l2, 1], maxcoeff=10**15, maxsteps=10**6)
    print(f"   remainder (P'''_arc - Z_e z3)/L_e on basis without z3: {rel2}")
    info['remainder_rel_without_z3'] = rel2
    out['arcs'].append(info)
print("total P'''(1) =", mp.nstr(tot3, 40))
out['P3_total'] = mp.nstr(tot3, dps)
json.dump(out, open(f'research/10_frontier/E1_edge_law/verify/out/t3_b{str(b).replace("/","_")}_{dps}.json', 'w'), indent=1, default=str)
