import sys, itertools
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *
h = Fr(1,2)
def cmul2(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def sadd2(*Ss):
    out = {}
    for S in Ss:
        for w, c in S.items():
            cur = out.get(w, (Fr(0),Fr(0)))
            out[w] = (cur[0]+c[0], cur[1]+c[1])
    return {w: c for w, c in out.items() if c != (Fr(0),Fr(0))}
def sscale2(S, cf):
    cf = (Fr(cf[0]), Fr(cf[1]))
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
e0 = (Fr(1),Fr(0)); e1 = (Fr(0),Fr(1))
b_vec = (Fr(0),Fr(-1)); om_b = (h,Fr(-1)); om_mb = (h,Fr(0))
r_vec = (-h,Fr(0)); om_r = (-h,Fr(-1))
def li_sym2(n, rv, omrv, cf):
    return sscale2({(lneg(omrv),) + (rv,)*(n-1): (Fr(1),Fr(0))}, cf)
L_vec = { (e1,): (Fr(1),Fr(0)) }
lam_vec = { (e0,): (Fr(1),Fr(0)) }
Li3_b  = li_sym2(3, b_vec, om_b, (1,0))
Li3_mb = li_sym2(3, b_vec, om_mb, (1,0))
Li3_r  = li_sym2(3, r_vec, om_r, (1,0))
Li2_b  = li_sym2(2, b_vec, om_b, (1,0))
P1 = sshuf2(L_vec, Li3_b)
P2 = sshuf2(lam_vec, Li3_r)
P3 = sshuf2(Li2_b, Li2_b)
L2 = sshuf2(L_vec, L_vec)
P4 = sshuf2(Li2_b, L2)
P5 = sshuf2(L2, L2)
for nm, P in [('L*Li3(b)', P1), ('lam*Li3(r)', P2), ('Li2(b)^2', P3), ('Li2(b)*L^2', P4), ('L^4', P5)]:
    d = delta22_spec(P, 2)
    print(f"{nm:>14}: delta22 = {d if d else 0}")
