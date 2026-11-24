import sys, json, re, itertools
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
h = Fr(1,2)
def cmul2(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def paircf(cf): return (Fr(cf[0]), Fr(cf[1]))

# symbols with (p,q) coefficients; words over letter-vectors in basis (e0=[2], e1=[1+sqrt2])
def sadd2(*Ss):
    out = {}
    for S in Ss:
        for w, c in S.items():
            cur = out.get(w, (Fr(0),Fr(0)))
            out[w] = (cur[0]+c[0], cur[1]+c[1])
    return {w: c for w, c in out.items() if c != (Fr(0),Fr(0))}
def sscale2(S, cf):
    return {w: cmul2(cf, c) for w, c in S.items() if cmul2(cf, c) != (Fr(0),Fr(0))}
def sshuf2(S, T):
    out = {}
    for u, cu in S.items():
        for v, cv in T.items():
            for w in shuffles(u, v):
                cc = cmul2(cu, cv)
                cur = out.get(w, (Fr(0),Fr(0)))
                out[w] = (cur[0]+cc[0], cur[1]+cc[1])
    return {w: c for w, c in out.items() if c != (Fr(0),Fr(0))}

# letters
e0 = (Fr(1),Fr(0)); e1 = (Fr(0),Fr(1))
b_vec = (Fr(0),Fr(-1)); om_b = (h,Fr(-1))        # [b] = -e1, [1-b] = e0/2 - e1
om_mb = (h,Fr(0))                                # [1+b] = e0/2
b2_vec = (Fr(0),Fr(-2)); om_b2 = (-h,Fr(1)); om_mb2 = (Fr(3,2),Fr(-1))  # [1-b^2]=? see below
r_vec = (-h,Fr(0)); om_r = (-h,Fr(-1))           # 1/sqrt2
half_vec = (Fr(-1),Fr(0)); om_half = (Fr(-1),Fr(0))
r2_vec = (-h,Fr(-1)); om_r2 = (-h,Fr(0))         # 1-1/sqrt2, 1-(1-1/sqrt2)=1/sqrt2
d_vec = (h,Fr(-1)); om_d = (Fr(0),Fr(-1))        # 2-sqrt2, 1-(2-sqrt2) = sqrt2-1 = b -> (0,-1)
# check 1-b^2: b^2 = 3-2sqrt2, 1-b^2 = 2sqrt2-2 = 2(√2-1) = 2/(1+sqrt2): [1-b^2] = (1,-1) -> -[1-b^2] = (-1,1) = om_b2 ✓

def li_sym2(n, rv, omrv, cf):
    return sscale2({(lneg(omrv),) + (rv,)*(n-1): (Fr(1),Fr(0))}, cf)

L_vec = { (e1,): (Fr(1),Fr(0)) }
lam_vec = { (e0,): (Fr(1),Fr(0)) }

# Li3 symbols
Li3_b  = li_sym2(3, b_vec, om_b, (1,0))
Li3_mb = li_sym2(3, b_vec, om_mb, (1,0))
Li3_r  = li_sym2(3, r_vec, om_r, (1,0))
Li2_b  = li_sym2(2, b_vec, om_b, (1,0))

# paper products
P1 = sshuf2(L_vec, sadd2(sscale2(Li3_b, (40,7)), sscale2(Li3_mb, (-40,7)), sscale2(Li3_r, (-72,0))))
P2 = sshuf2(lam_vec, sadd2(sscale2(Li3_b, (6,13)), sscale2(Li3_mb, (6,-13)), sscale2(Li3_r, (Fr(-47,2)*0,0))))
# careful: -47*sqrt2/2 * Li3(1/sqrt2): coefficient (0, -47/2)
P2 = sshuf2(lam_vec, sadd2(sscale2(Li3_b, (6,13)), sscale2(Li3_mb, (6,-13)), sscale2(Li3_r, (0,-Fr(47,2)))))
P3 = sscale2(sshuf2(Li2_b, Li2_b), (30,0))
L2 = sshuf2(L_vec, L_vec)
lamL = sshuf2(lam_vec, L_vec)
lam2 = sshuf2(lam_vec, lam_vec)
P4 = sshuf2(Li2_b, sadd2(sscale2(L2, (30, Fr(7,2))), sscale2(lamL, (-24,0)), sscale2(lam2, (0, -Fr(19,8)))))
P5 = sadd2(sscale2(sshuf2(L2, L2), (Fr(13,24), -Fr(25,24))),   # L^4
           sscale2(sshuf2(sshuf2(L2, L_vec), lam_vec), (-Fr(73,6), -Fr(29,3))),  # L^3 lam
           sscale2(sshuf2(L2, lam2), (Fr(63,8), Fr(67,16))),   # L^2 lam^2
           sscale2(sshuf2(L_vec, sshuf2(lam2, lam_vec)), (Fr(41,24), Fr(17,6))), # L lam^3
           sscale2(sshuf2(lam2, lam2), (Fr(-61,192), Fr(137,192))))  # lam^4
SW4prods = sadd2(P1, P2, P3, P4, P5)
print("paper product symbol terms:", len(SW4prods))
dP = delta22_spec(SW4prods, 2)
print("delta22(paper products):", dP if dP else "= 0 (must be)")
