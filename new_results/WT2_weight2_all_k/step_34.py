code = r'''
"""Build M_k integrand PF in v from the closed form (genintegrand + fastred)."""
from fractions import Fraction as Fr
from math import factorial
from genintegrand import S_poly, p_comp_1mt, pmul as qpmul, ptrim
from fastred import *

SQ2Q = QSI(0,1,0,0)   # sqrt2

def qpoly_to_qsi(P):
    return [QSI(c) for c in P]

def qsi_pmul(P, Q):
    out = [ZERO]*(len(P)+len(Q)-1)
    for a,ca in enumerate(P):
        for b,cb in enumerate(Q):
            out[a+b] = out[a+b] + ca*cb
    while out and out[-1].is_zero(): out.pop()
    return out

def qsi_pscale(P, q):
    return [c*q for c in P]

TV = [QSI(-1), QSI(2), QSI(1)]    # v^2 + 2v - 1
FOURV = [ZERO, QSI(4)]            # 4v
V2P1 = [QSI(1), ZERO, QSI(1)]     # v^2+1
VMINUSA = [-AQ, ONE]
VPLUSA = [AQ, ONE]
VMINUSB = [-BQ, ONE]
VPLUSB = [BQ, ONE]
VPOLY = [ZERO, ONE]

def Stilde(S):
    """S(t(v)) * (4v)^D as polynomial in v."""
    D = len(S)-1
    out = [ZERO]
    for s, c in enumerate(S):
        if c == 0: continue
        term = [ONE]
        for _ in range(s):
            term = qsi_pmul(term, TV)
        for _ in range(D-s):
            term = qsi_pmul(term, FOURV)
        term = qsi_pscale(term, QSI(c))
        # add
        n = max(len(out), len(term))
        out += [ZERO]*(n-len(out)); term += [ZERO]*(n-len(term))
        out = [a+b for a,b in zip(out, term)]
    while out and out[-1].is_zero(): out.pop()
    return out, D

def sqrt_part_PF(k, S):
    """PF of  -sqrt2 8^k 4^{-D} v^{2k-1-D} Stilde(v) (v^2+1)^{1-2k} / ((v^2-a^2)(v^2-b^2))."""
    St, D = Stilde(S)
    if not St:   # S == 0
        return pf_zero()
    e = 2*k - 1 - D
    num = qsi_pscale(St, SQ2Q * QSI(Fr(-(8**k), 4**D)))
    den_extra_v = 0
    if e >= 0:
        for _ in range(e):
            num = qsi_pmul(num, VPOLY)
    else:
        den_extra_v = -e
    den = qsi_pmul(qsi_pmul(VMINUSA, VPLUSA), qsi_pmul(VMINUSB, VPLUSB))
    for _ in range(den_extra_v):
        den = qsi_pmul(den, VPOLY)
    for _ in range(2*k-1):
        den = qsi_pmul(den, V2P1)
    return pf_from_poly(num, den)

def AB_part_PF(k):
    """-(v^2+1)/(2 k! v (v-b)(v+a)) for each of A and folded B."""
    den = qsi_pmul(qsi_pmul(VPOLY, VMINUSB), VPLUSA)
    num = qsi_pscale(V2P1, QSI(Fr(-1, 2*factorial(k))))
    return pf_from_poly(num, den)

def build_integrand_pf(k):
    """{key (i,j), i>=j, i+j=k: PF} for M_k = int_b^a sum R_ij l1^i l2^j dv."""
    S = {}
    for m in range(0, k+1):
        S[m] = S_poly(k, m)
    out = {}
    for i in range(k, (k-1)//2, -1):
        j = k - i
        P = sqrt_part_PF(k, S[i])            # unfolded l1^i l2^j
        Sj = p_comp_1mt(S[j])                # folded l1^j l2^i  (t -> 1-t)
        P = pf_add(P, sqrt_part_PF(k, Sj))
        if j == 0:                           # key (k,0): A-part + folded B-part
            P = pf_add(P, AB_part_PF(k))
            P = pf_add(P, AB_part_PF(k))
        out[(i, j)] = P
    return out
'''
with open('buildintegrand.py','w') as f:
    f.write(code)
print("written")
