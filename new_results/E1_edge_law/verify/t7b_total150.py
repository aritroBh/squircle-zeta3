# T7b: PSLQ on the TOTAL P'''(1) of the rational hexagon (b=0), second precision (150 digits)
import json, mpmath as mp
mp.mp.dps = 150
d = json.load(open('research/10_frontier/E1_edge_law/verify/out/t3_b0_150.json'))
Pt = mp.mpf(d['P3_total'])
z3, l2, pi2 = mp.zeta(3), mp.log(2), mp.pi**2
s2, s5 = mp.sqrt(2), mp.sqrt(5)
B = [z3, pi2*l2, l2**3, pi2, l2**2, l2, mp.mpf(1)]
rel = mp.pslq([Pt] + [s5*b for b in B] + [s2*b for b in B], maxcoeff=10**12, maxsteps=10**7)
print("TOTAL P'''(1), 150 digits, PSLQ on [P''', sqrt5*(z3,pi2 log2,log2^3,pi2,log2^2,log2,1), sqrt2*(...)]:", rel)
print(f"total zeta(3) coefficient = ({-rel[1]}/{rel[0]}) sqrt5 + ({-rel[8]}/{rel[0]}) sqrt2 ; edge law = (72/125) sqrt5 - (15/64) sqrt2")
print("residual:", mp.nstr(mp.fsum(r*x for r, x in zip(rel, [Pt] + [s5*b for b in B] + [s2*b for b in B])), 3))
Rt = Pt - (72*s5/125 - 15*s2/64)*z3
rel2 = mp.pslq([Rt] + [s5*b for b in B[1:]] + [s2*b for b in B[1:]], maxcoeff=10**12, maxsteps=10**7)
print("P''' - (sum Z_e) z3 on z3-free basis {sqrt5,sqrt2}x6 (maxcoeff 1e12):", rel2)
