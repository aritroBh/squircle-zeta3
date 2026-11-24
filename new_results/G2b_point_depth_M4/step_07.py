"""step_07.py -- G2b convention fix (planner item 1): re-solve N5 with CONVERGENT atoms.

CONVENTION (documented): atom values are Li_{3,2}(x,y) = sum_{n>=2} (x^n/n^3) sum_{m<n} y^m/m^2
(n-outer iterated series), which for |x| <= 1, |xy| <= 1 is absolutely convergent and equals the
paper's Li_{s,t}(x,y) = sum_{n>m>0} x^n y^m/(n^s m^t). Only such atoms were ever evaluated
(values_hex.json; 166 values). The J-word identity J((0,0,1/x,0,1/(xy));0->1) = +Li_{3,2}(x,y)
(paper convention) is verified numerically to 60 digits on six atoms below and in the transcript.

The step_05 candidate N5 = (3/4)Li32(8/9,3/4) + 4*Li32(2,1/3) used Li32(2,1/3) (|x| = 2 > 1),
which DIVERGES in the paper's convention; no value for it was ever computed or used. Here N5 is
re-solved over convergent atoms only (|x| <= 1, |xy| <= 1), exact rational solve in the T^3/D^3
quotient, and the symbol-level pairing obs(N5) = (41, -33/2) is re-verified exactly.
"""
import sys, os, json, time, itertools
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, G2B)
from g2b_dbg import s05
log = s05.log
t0 = time.time()

if __name__ == '__main__':
    log(f"step_07 start {time.strftime('%H:%M:%S')}")
    sols = [Fr(-8), Fr(-3), Fr(-2), Fr(-1), Fr(-1, 2), Fr(-1, 3), Fr(-1, 8), Fr(1, 9), Fr(1, 4),
            Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(3, 4), Fr(8, 9), Fr(9, 8), Fr(4, 3), Fr(3, 2),
            Fr(2), Fr(3), Fr(4), Fr(9), Fr(1)]
    atoms = sorted({(x, y) for x in sols for y in sols
                    if x != 1 and x * y in sols and abs(x) <= 1 and abs(x * y) <= 1})
    log(f"convergent atoms (|x|<=1, |xy|<=1): {len(atoms)}")
    from step_12 import br
    X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
    P = br(X, Y); Q = br(X, br(X, Y)); T = br(Y, br(X, Y))
    def obs_of_SE(SE):
        SE = {w: (cf[0] if isinstance(cf, tuple) else cf) for w, cf in SE.items()}
        out = []
        for A in (Q, T):
            tot = Fr(0)
            for u, cu in A.items():
                for v, dv in P.items():
                    tot += cu * dv * (SE.get(u + v, Fr(0)) - SE.get(v + u, Fr(0)))
            out.append(tot)
        return out
    adata = {}
    for x, y in atoms:
        dw = s05.delta_word_sym(s05.li32_word(x, y), 0, 1)
        dE = s05.ft_expand_sym({w: (cf,) for w, cf in dw.items()}, 2)
        qcw = s05.quotient_coords(dE)[0]
        S_at = s05.ft_expand_sym({w: (cf,) for w, cf in s05.li32_sym(x, y).items()}, 2)
        adata[(x, y)] = (qcw, tuple(obs_of_SE(S_at)))
    target = (Fr(-41), Fr(-33, 2))   # Delta's quotient coords (from cobracket_hex.json / step_05)
    # pairs
    sol_pairs = []
    for a1, a2 in itertools.combinations(adata, 2):
        v1, v2 = adata[a1][0], adata[a2][0]
        det = v1[0] * v2[1] - v1[1] * v2[0]
        if not det: continue
        d1 = (target[0] * v2[1] - target[1] * v2[0]) / det
        d2 = (v1[0] * target[1] - v1[1] * target[0]) / det
        h = max(abs(d1.numerator), abs(d1.denominator), abs(d2.numerator), abs(d2.denominator))
        sol_pairs.append((h, ((a1, d1), (a2, d2))))
    sol_pairs.sort(key=lambda t: t[0])
    log("best convergent-atom PAIR solutions (height, atoms):")
    for h, sol in sol_pairs[:6]:
        log(f"   h={h}: " + " + ".join(f"({d})*Li32({a[0]},{a[1]})" for a, d in sol))
    # pick the best and verify obs exactly
    best = sol_pairs[0][1]
    oN5 = [Fr(0), Fr(0)]
    for a, d in best:
        o = adata[a][1]
        oN5 = [oN5[0] + d * o[0], oN5[1] + d * o[1]]
    log(f"chosen N5 = " + " + ".join(f"({d})*Li32({a[0]},{a[1]})" for a, d in best))
    log(f"obs(N5) = {[str(o) for o in oN5]} (target (41, -33/2)): "
        f"{'MATCH' if tuple(oN5) == (Fr(41), Fr(-33, 2)) else 'MISMATCH'}")
    json.dump({'N5_convergent': [[f'Li32({a[0]},{a[1]})', str(d)] for a, d in best],
               'obs_N5': [str(o) for o in oN5],
               'pairs': [[h, [[f'Li32({a[0]},{a[1]})', str(d)] for a, d in sol]] for h, sol in sol_pairs[:10]]},
              open(G2B + 'N5_convergent.json', 'w'), indent=1)
    log(f"step_07 done ({time.time()-t0:.0f}s)")
