# T6: PSLQ directly on the finite-difference values (no d-expansion anywhere), at two precisions.
import json, mpmath as mp
runs = [('t4_fd_120_0.001_20.json', 38), ('t4_fd_170_0.0004_26.json', 60)]
for fn, digs in runs:
    mp.mp.dps = digs
    d = json.load(open('research/10_frontier/E1_edge_law/verify/out/' + fn))
    z3, l2, pi2 = mp.zeta(3), mp.log(2), mp.pi**2
    s5, s2 = mp.sqrt(5), mp.sqrt(2)
    P3 = mp.mpf(d['arc1'][3])
    X = P3/s5
    rel = mp.pslq([X, z3, pi2*l2, l2**3, pi2, l2**2, l2, 1], maxcoeff=10**6, maxsteps=10**6)
    print(f"== {fn}, working at {digs} digits")
    print(f"   arc1 (r=2 edge): PSLQ [P'''/sqrt5, z3, pi2 log2, log2^3, pi2, log2^2, log2, 1] -> {rel}")
    if rel: print(f"   zeta(3) coefficient of P'''_arc1 = ({-rel[1]}/{rel[0]}) sqrt5 ; edge law Z_e = (72/500) sqrt5")
    Z = 18*s5/125
    R = (P3 - Z*z3)/s5
    rel2 = mp.pslq([R, pi2*l2, l2**3, pi2, l2**2, l2, 1], maxcoeff=10**6, maxsteps=10**6)
    print(f"   remainder (P'''_arc1 - Z_e z3)/sqrt5 on the z3-free basis (maxcoeff 1e6): {rel2}")
    # the total perimeter: P''' = 4*arc1 + 2*arc3 ; test it against the edge law with the full 14-element basis
    Pt = mp.mpf(d['total'][3])
    Zt = 72*s5/125 - 15*s2/64
    Rt = Pt - Zt*z3
    B = [pi2*l2, l2**3, pi2, l2**2, l2, 1]
    rel3 = mp.pslq([Rt] + [s5*b for b in B] + [s2*b for b in B], maxcoeff=10**6, maxsteps=10**6)
    print(f"   TOTAL perimeter: (P''' - sum_e Z_e z3) on {{sqrt5,sqrt2}} x z3-free basis: {rel3}")
    rel4 = mp.pslq([Pt] + [s5*z3, s2*z3] + [s5*b for b in B] + [s2*b for b in B], maxcoeff=10**6, maxsteps=10**6)
    print(f"   TOTAL perimeter: P''' on {{sqrt5,sqrt2}} x full basis (z3 first): {rel4}")
    if rel4: print(f"   => total zeta(3) coefficient = ({-rel4[1]}/{rel4[0]}) sqrt5 + ({-rel4[2]}/{rel4[0]}) sqrt2 ; edge law = (72/125) sqrt5 - (15/64) sqrt2")
