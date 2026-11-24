# s14: TRANSLATION check for the adjacent identities.  prove_adj56.feed_terms(delta,s) is the claimed Newton-space RHS
# (A-parent, diagonal direct DL2, fold).  Here the SAME three parts are computed from T1 itself at integer p (true binomials,
# actual diagonal tower, K=k-2p=1..5), extrapolated to a polynomial in K and converted to the Newton basis binom(K,s);
# the script's terms (evaluated at integer p) must agree PART BY PART.
import sys, time
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from math import comb
import sympy as sp
from candlib import *
import prove_adj56 as PA
from t1lib import _mi_pow
t0 = time.time()
cand = Cand(base_exprs())
taus = cand.taus
def parts_at(p, K, delta):
    q = p-1; d = p+q; k = 2*p+K
    A = ZERO
    for dp in range(delta+1):
        A = A + cand.get(p+1, q, k, dp).mul(taus.t1[delta-dp] * Fr(p+1, 2*d-dp-1))
    D = ZERO
    for dp in range(delta+1):
        D = D + cand.get(p, p, k, dp).mul(taus.t2[delta-dp] * Fr(q+1, 2*d-dp-1))
    mm = 4*p-4-delta
    Fo = ZERO
    for mu in range(mm, 4*p-4+1):
        g = ZERO
        for dp in range(0, 4*p-4-mu+1):
            g = g - cand.get(p, p, k, dp).mul(taus.t1[4*p-4-mu-dp] / GQ(4*p-3-dp))
        f = _mi_pow(mu) * ipow(mm-2) * comb(mu-2, mu-mm)
        Fo = Fo + g.conj().mul(f * (-p))
    e0 = E00(p, q)
    out = []
    for X in (A, D, Fo):
        assert X.x.iszero()
        out.append(X.y / e0.y)
    return out        # GQ ratios to rho_{p,p-1}(0)
def lagrange_extrap(vals, xs, x):
    tot = GQ(0)
    for i_, xi in enumerate(xs):
        num = Fr(1); den = Fr(1)
        for j_, xj in enumerate(xs):
            if i_ == j_: continue
            num *= (x - xj); den *= (xi - xj)
        tot = tot + vals[i_] * (num/den)
    return tot
def to_gq(z):
    z = sp.nsimplify(sp.expand(z)); re_, im_ = z.as_real_imag()
    return GQ(Fr(int(re_.p), int(re_.q)), Fr(int(im_.p), int(im_.q)))
pS = sp.Symbol('p')
ntest = 0; nbad = 0
for p in (3, 4, 5, 6, 7, 9, 12):
    for delta in (5, 6):
        P = [parts_at(p, K, delta) for K in range(1, 6)]            # K = 1..5, three parts each
        for part in range(3):
            vals = [P[K-1][part] for K in range(1, 6)]
            # polynomial of degree <= floor(delta/2) <= 3 through K=1..4; K=5 is a degree check
            chk = lagrange_extrap(vals[:4], [1, 2, 3, 4], 5)
            assert chk == vals[4], ('degree>3', p, delta, part)
            ext = [lagrange_extrap(vals[:4], [1, 2, 3, 4], j) for j in range(0, 4)]
            for s in range(0, delta//2+1):
                newton = GQ(0)
                for j in range(0, s+1):
                    newton = newton + ext[j] * ((-1)**(s-j) * comb(s, j))
                terms = PA.feed_terms(delta, s)
                sl = terms[part*(delta+1):(part+1)*(delta+1)]
                script = to_gq(sum(t.subs(pS, p) for t in sl))
                ntest += 1
                if script != newton:
                    nbad += 1; print('MISMATCH', p, delta, part, s, script, newton)
    print(f'  p={p} done ({round(time.time()-t0,1)}s)', flush=True)
print(f'part-by-part comparison (A-parent / diagonal-direct / fold) of prove_adj56.feed_terms vs T1-derived Newton coefficients: {ntest} comparisons, {nbad} mismatches')
