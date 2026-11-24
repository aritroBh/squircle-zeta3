"""check_T1_large.py -- referee check that the closed-form e/ehat families satisfy the ENGINE's T1 operator
(TOWN1 trunc_engine.Lambda_ij, tau constants computed from QSI definitions) at degrees FAR beyond the data
(d up to ~500), i.e. that tow3lib.Lam / prove_adj56.feed_terms are faithful transcriptions of T1 and the symbolic
identities really are statements about the tower recursion.  Child normalised to e00 = 1; parents scaled by
the e00 ratios used in the proofs (verified against data for p<=26 in check_data.py)."""
import sys, random, time
from fractions import Fraction as Fr
from math import comb
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/08_theory/TOWN1_tower_large_n')
import sympy as sp
import trunc_engine as te
from fastred import QSI
import check_data as CD
import tow3lib as T3
import prove_E3b56 as PE

def G(re, im=0): return QSI(Fr(re), 0, Fr(im), 0)
def gq(pair): return G(pair[0], pair[1])

evE = {}; evH = {}
for dl in range(0, 7):
    for s in range(0, dl//2+1):
        evE[(dl, s)] = CD.make_eval(T3.e(dl, s), (T3.d, T3.m))
        evH[(dl, s)] = CD.make_eval(sp.cancel(sp.together(PE.ehat(dl, s))), (PE.p,))

def e_val(dl, s, dd, mm):
    if s > dl//2: return G(0)
    return gq(evE[(dl, s)](dd, mm))
def eh_val(dl, s, pp):
    if s > dl//2: return G(0)
    return gq(evH[(dl, s)](pp))

def poly_off(dl, dd, mm, k):      # sum_s C(k-dd, s) e(dl,s)(dd,mm), off-diagonal key of degree dd
    n = k - dd; assert n >= 0
    tot = G(0)
    for s in range(dl//2+1): tot = tot + e_val(dl, s, dd, mm).scale(comb(n, s))
    return tot
def poly_diag(dl, pp, k):          # stable diagonal (pp,pp): normalised by rho_pp(0)
    if k == 2*pp:
        return e_val(dl, 0, 2*pp, pp)            # initial level = rho_pp(dl)/rho_pp(0) = a_dl(2pp,pp)
    n = k - 2*pp - 1; assert n >= 0
    tot = G(0)
    for s in range(dl//2+1): tot = tot + eh_val(dl, s, pp).scale(comb(n, s))
    return tot

I_ = G(0, 1)
def run():
    random.seed(20261028)
    t0 = time.time()
    DL = 6
    tot = {'nonadj': 0, 'adj': 0, 'diag': 0}; bad = {'nonadj': 0, 'adj': 0, 'diag': 0}
    # (a) non-adjacent children
    keys = []
    for _ in range(25):
        q = random.choice([0, 0, 1, 2, random.randint(3, 150)])
        p = q + random.randint(2, 300)
        keys.append((p, q))
    keys += [(2+0, 0), (5, 0), (5, 1), (6, 2), (7, 0), (9, 3)]
    for (p, q) in keys:
        d = p+q
        if 2*d-4 < DL: continue
        rA = I_.scale(Fr(2*d-1, p+1)); rB = I_.scale(Fr(-(2*d-1), q+1))
        for k in range(d+1, d+10):
            def U(key, dp, _k=k, _p=p, _q=q, _d=d, _rA=rA, _rB=rB):
                if key == (_p+1, _q): return _rA * poly_off(dp, _d+1, _p+1, _k)
                if key == (_p, _q+1): return _rB * poly_off(dp, _d+1, _p, _k)
                raise KeyError(key)
            for dl in range(0, DL+1):
                got = te.Lambda_ij(p, q, dl, U)
                pred = poly_off(dl, d, p, k)
                tot['nonadj'] += 1
                if not (got - pred).is_zero():
                    bad['nonadj'] += 1
                    if bad['nonadj'] <= 3: print('NONADJ BAD', p, q, dl, k)
    # (b) adjacent children (p,p-1), diagonal parent via ehat; includes the fold
    ps = sorted(set([3, 4, 5, 6, 7, 10, 25, 26, 27, 40, 77, 123, 250] + [random.randint(28, 400) for _ in range(8)]))
    for p in ps:
        q = p-1; d = 2*p-1
        if 2*d-4 < DL: continue
        rA = I_.scale(Fr(2*d-1, p+1))
        rD = G(0, Fr(-(4*p-3), 2*p))        # rho_pp(0)/e00_child
        for k in range(d+1, d+10):
            def U(key, dp, _k=k, _p=p, _q=q, _d=d, _rA=rA, _rD=rD):
                if key == (_p+1, _q): return _rA * poly_off(dp, _d+1, _p+1, _k)
                if key == (_p, _p): return _rD * poly_diag(dp, _p, _k)
                raise KeyError(key)
            for dl in range(0, DL+1):
                got = te.Lambda_ij(p, q, dl, U)
                pred = poly_off(dl, d, p, k)
                tot['adj'] += 1
                if not (got - pred).is_zero():
                    bad['adj'] += 1
                    if bad['adj'] <= 3: print('ADJ BAD', p, q, dl, k)
    # (c) diagonal children (p,p), parent (p+1,p)
    for p in sorted(set(ps + [2, 3, 4, 5])):
        d = 2*p
        if 2*d-4 < 0: continue
        rP = I_.scale(Fr(2*(4*p-1), p+1)) if False else G(0, Fr(2*(4*p-1), p+1))   # e00_{(p+1,p)}/rho_pp(0)
        for k in range(d+1, d+10):
            def U(key, dp, _k=k, _p=p, _rP=rP):
                if key == (_p+1, _p): return _rP * poly_off(dp, 2*_p+1, _p+1, _k)
                raise KeyError(key)
            for dl in range(0, min(DL, 2*d-4)+1):
                got = te.Lambda_ij(p, p, dl, U)
                pred = poly_diag(dl, p, k)
                tot['diag'] += 1
                if not (got - pred).is_zero():
                    bad['diag'] += 1
                    if bad['diag'] <= 3: print('DIAG BAD', p, dl, k)
    print('engine-T1 vs closed-form families:', tot, 'bad:', bad, 'time', round(time.time()-t0, 1))
run()
