
"""Degree-d kernels from the degree-(d+1) tower (WTS Lemma 2 / TOW1), exact QSI."""
import os, sys
from fractions import Fraction as Fr
_HERE = os.path.dirname(os.path.abspath(__file__))
_WTS = os.path.normpath(os.path.join(_HERE, '..', 'WTS_weight_k_structure', 'theory'))
sys.path.insert(0, _WTS)
from fastred import (QSI, ZERO, ONE, IQ, AQ, BQ, pf_zero, pf_add, pf_scale,
                     pf_integrate_higher, pf_eval_at, pf_mul_DL, pf_split, pf_fold, DL1, DL2)

def unq4(t): return QSI(*t)

def gm_from_tower(h_orders):
    """h_orders: dict {r: QSI} coefficient of (v-i)^{-r}, r>=2 (the +i higher part).
       Returns Gm = (1 - eval_b) o Int(HigherPart) as a PF [real, ipos, ineg, const]."""
    maxr = max(h_orders.keys())
    ipos = [ZERO]  # order 1 = 0
    ineg = [ZERO]
    for r in range(2, maxr+1):
        c = h_orders.get(r, ZERO)
        ipos.append(c)
        ineg.append(c.conj_i())
    H = pf_zero(); H[1] = ipos; H[2] = ineg
    G = pf_integrate_higher(H)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], ZERO - C]

def simple_part_of(W):
    simple, higher, const = pf_split(W)
    assert const.is_zero()
    return simple, higher

def kernel_20_21_level(h30, h21):
    """degree-2 kernels at a level, from degree-3 towers h30, h21 (dicts {r: QSI}).
       Returns (K20_simple, K11_simple, W20_higher, W11_higher)."""
    gm30 = gm_from_tower(h30)
    gm21 = gm_from_tower(h21)
    W20 = pf_add(pf_scale(pf_mul_DL(gm30, DL1), Fr(-3)),
                 pf_scale(pf_mul_DL(gm21, DL2), Fr(-1)))
    W11 = pf_scale(pf_mul_DL(gm21, DL1), Fr(-2))
    K20, Hh20 = simple_part_of(W20)
    K11, Hh11 = simple_part_of(W11)
    return K20, K11, Hh20, Hh11

def kernel_deg3_level(h40, h31, h22):
    """degree-3 kernels K30, K21 from degree-4 towers (dicts {r: QSI})."""
    gm40 = gm_from_tower(h40)
    gm31 = gm_from_tower(h31)
    gm22 = gm_from_tower(h22)
    W30 = pf_add(pf_scale(pf_mul_DL(gm40, DL1), Fr(-4)),
                 pf_scale(pf_mul_DL(gm31, DL2), Fr(-1)))
    # W21 = -3 Gm31 DL1 - 2 Gm22 DL2 - 2 F(Gm22 DL1)
    t1 = pf_scale(pf_mul_DL(gm31, DL1), Fr(-3))
    t2 = pf_scale(pf_mul_DL(gm22, DL2), Fr(-2))
    t3 = pf_scale(pf_fold(pf_mul_DL(gm22, DL1)), Fr(-2))
    W21 = pf_add(pf_add(t1, t2), t3)
    K30, Hh30 = simple_part_of(W30)
    K21, Hh21 = simple_part_of(W21)
    return K30, K21, Hh30, Hh21

# ---- closed forms (WT3) ----
SQ2 = QSI(0,1,0,0)
def pf_real(terms):
    """terms: list of (pole QSI real, coeff QSI). Returns PF simple (real dict only)."""
    P = pf_zero()
    for r, c in terms:
        P[0][r] = P[0].get(r, ZERO) + c
    for r in list(P[0].keys()):
        if P[0][r].is_zero(): del P[0][r]
    return P

def R20_closed(k):
    # (1-sqrt2)/(4v) + 1/(4(v-a)) + (2k-7)/(4(v+a)) - (k-2)/(2(v+b))
    return pf_real([
        (ZERO, (ONE - SQ2).scale(Fr(1,4))),
        (AQ,   ONE.scale(Fr(1,4))),
        (-AQ,  QSI(2*k-7).scale(Fr(1,4))),
        (-BQ,  QSI(k-2).scale(Fr(-1,2))),
    ])

def R11_closed(k):
    # (2k-5+sqrt2)/(4v) - (2k-5)/(2(v+a))
    return pf_real([
        (ZERO, (QSI(2*k-5) + SQ2).scale(Fr(1,4))),
        (-AQ,  QSI(2*k-5).scale(Fr(-1,2))),
    ])

def pf_equal(P, Q):
    # compare real dicts and ipos/ineg/const
    for r in set(P[0]) | set(Q[0]):
        if P[0].get(r, ZERO) != Q[0].get(r, ZERO): return False
    def lst(L):
        while L and L[-1].is_zero(): L.pop()
        return L
    if lst(list(P[1])) != lst(list(Q[1])): return False
    if lst(list(P[2])) != lst(list(Q[2])): return False
    return P[3] == Q[3]
