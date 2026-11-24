"""Drinfeld associator Phi(e0,e1) (BSST arXiv:1304.7304 conventions) at 5 points, first column,
as Q-vectors over the engine's MZV basis BASIS[W] (reduction via the repo's read-only DSTable).
Convention (BSST eq.5,10): zeta_{n1..nr} = int_{0<z1<..<zn<1} w1 w0^{n1-1} ... w1 w0^{nr-1}  (SS ordering k1<..<kr)
   Phi = sum_w  rev(w)[e0,e1] * zeta^{shuffle-reg}(w),   zeta(w0)=zeta(w1)=0.
In the repo (mzv_exact) convention zeta(s1..sd)=sum n1>..>nd, so SS (n1..nr) == repo (nr..n1)."""
import sys, os
sys.dont_write_bytecode = True
CHN = "<project>/research/08_theory/CHN1_chain_bodies"
sys.path.insert(0, CHN)
from fractions import Fraction as F
from functools import lru_cache
import itertools
import mzv_exact as M
import engine as E

# ---------------- shuffle regularisation of words over {0,1} ----------------
def insert_all(w, ch):
    """all (position, word) with ch inserted into w at position p=0..len(w)"""
    return [(p, w[:p] + ch + w[p:]) for p in range(len(w) + 1)]

@lru_cache(maxsize=None)
def reg(w):
    """zeta^{sh}(w) = sum coeff * zeta(convergent word); convergent: starts '1', ends '0' (or empty)."""
    if w == '':
        return {'': F(1)}
    if w[0] == '0':
        a = len(w) - len(w.lstrip('0'))
        wp = w[1:]
        out = {}
        for p, ww in insert_all(wp, '0'):
            if p >= a:                      # positions p<a reproduce w itself (a copies)
                for k, c in reg(ww).items():
                    out[k] = out.get(k, F(0)) - c / a
        return {k: c for k, c in out.items() if c}
    if w[-1] == '1':
        b = len(w) - len(w.rstrip('1'))
        wp = w[:-1]
        out = {}
        for p, ww in insert_all(wp, '1'):
            if p < len(wp) - (b - 1):       # positions in trailing-ones block reproduce w (b copies)
                for k, c in reg(ww).items():
                    out[k] = out.get(k, F(0)) - c / b
        return {k: c for k, c in out.items() if c}
    return {w: F(1)}

def conv_to_repo_index(cw):
    """convergent word '1 0^{n1-1} 1 0^{n2-1} ...' -> repo-convention tuple (nr,...,n1)."""
    ns = []
    i = 0
    while i < len(cw):
        assert cw[i] == '1'
        j = i + 1
        while j < len(cw) and cw[j] == '0':
            j += 1
        ns.append(j - i)
        i = j
    return tuple(reversed(ns))

def zeta_vec(w):
    """zeta^{sh}(w) as list over BASIS[|w|] (weight |w|>=2) -- Fractions."""
    W = len(w)
    T = E.table(W)
    acc = [F(0)] * T.nullity
    for cw, c in reg(w).items():
        idx = conv_to_repo_index(cw)
        assert idx[0] >= 2 and sum(idx) == W
        rv = T.reduce_word(idx)
        for i in range(T.nullity):
            acc[i] += c * rv[i]
    return acc

# ---------------- BSST five-point 6x6 matrices (eq. below (23)) ----------------
def e_mats(s12, s13, s23, s24, s34, s123, s234):
    z = 0
    e0 = [[s123, 0, -s13 - s23, -s12, -s12, s12],
          [0, s123, -s13, -s12 - s23, s13, -s13],
          [0, 0, s12, 0, -s12, 0],
          [0, 0, 0, s13, 0, -s13],
          [0] * 6, [0] * 6]
    e1 = [[0] * 6, [0] * 6,
          [s34, 0, -s34, 0, 0, 0],
          [0, s24, 0, -s24, 0, 0],
          [s34, -s34, s23 + s24, s34, -s234, 0],
          [-s24, s24, s24, s23 + s34, 0, -s234]]
    return e0, e1

def matvec(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]

def phi_first_column(e0, e1, Wmax, ring=F):
    """returns dict W -> list(6) of vectors over BASIS[W] : (Phi e_1)_i at weight W (W>=2); weight 0,1 handled separately."""
    e = {'0': e0, '1': e1}
    e1v = [ring(1)] + [ring(0)] * 5
    vecs = {'': e1v}
    out = {}
    cur = {'': e1v}
    for W in range(1, Wmax + 1):
        nxt = {}
        for w, v in cur.items():
            for a in '01':
                nxt[w + a] = matvec(e[a], v)     # v_{wa} = E_a v_w  (rev(w a) = a rev(w))
        cur = nxt
        if W >= 2:
            tot = [[F(0)] * E.table(W).nullity for _ in range(6)]
            for w, v in cur.items():
                if all(x == 0 for x in v):
                    continue
                zv = zeta_vec(w)
                for i in range(6):
                    if v[i] != 0:
                        for k in range(len(zv)):
                            tot[i][k] += v[i] * zv[k]
            out[W] = tot
    return out

# kinematic dictionary: our integral <-> F^{(23)} at s23=0
def kin_from_a(a):
    a1, a2, a3, a4 = [F(x) for x in a]
    s12 = a3; s23 = F(0); s34 = a2; s45 = a1 + a3; s51 = a2 + a4
    s13 = s45 - s12 - s23          # = a1
    s24 = s51 - s23 - s34          # = a4
    return dict(s12=s12, s13=s13, s23=s23, s24=s24, s34=s34, s123=s45, s234=s51)

if __name__ == '__main__':
    # self test: weight 2,3 structure of Phi
    print(reg('010'), reg('110'), reg('100'))
