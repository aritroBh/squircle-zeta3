# s09: FIRST-PRINCIPLES check of a_delta = rho(delta)/rho(0), delta<=6, from the TOW1 s.3 definition of S_m(t) and
# A(v) = v^(2d-1) S_m(t(v)) (v+i)^(2-2d) / (v^4-6v^2+1)  (exact truncated Taylor series at v=i over Q(i)); independent of
# the R_l lemmas, the Bell assembly and the TOWN1 data.  Compared with tow3lib.arho (the closed forms a_5, a_6, ...).
import sys, time
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp
import tow3lib as T3
from t1lib import GQ, RF
t0 = time.time()
N = 7
def smul(a, b):
    out = [GQ(0) for _ in range(N)]
    for i, x in enumerate(a):
        if x.iszero(): continue
        for j, y in enumerate(b):
            if i+j >= N: break
            out[i+j] = out[i+j] + x*y
    return out
def sadd(a, b): return [x+y for x, y in zip(a, b)]
def sscale(a, c): return [x*c for x in a]
def sconst(c): return [c if isinstance(c, GQ) else GQ(c)] + [GQ(0) for _ in range(N-1)]
def sinv(a):
    out = [GQ(0) for _ in range(N)]
    out[0] = GQ(1)/a[0]
    for n in range(1, N):
        acc = GQ(0)
        for k in range(1, n+1): acc = acc + a[k]*out[n-k]
        out[n] = -acc*out[0]
    return out
def spow(a, n):
    if n < 0: return spow(sinv(a), -n)
    r = sconst(1); b = a
    while n:
        if n & 1: r = smul(r, b)
        b = smul(b, b); n >>= 1
    return r
def stirling2(n, k):
    return int(sp.functions.combinatorial.numbers.stirling(n, k))
def w(s, m):
    return Fr(factorial(s)*stirling2(m, s)*2**m, factorial(m))
def binom_half(j):
    r = Fr(1)
    for u in range(j): r *= (Fr(1, 2)-u)
    return r/factorial(j)
I_ = GQ(0, 1)
v = [I_, GQ(1)] + [GQ(0)]*(N-2)
vinv = sinv(v)
t = sadd(sadd(sscale(v, Fr(1, 4)), sconst(Fr(1, 2))), sscale(vinv, Fr(-1, 4)))
omt = sadd(sconst(1), sscale(t, -1))
q2 = sadd(smul(t, t), smul(omt, omt))
def rho_ratios(dd, mm):
    # S_m(t(v)) as a series in s = v - i
    tp = {0: sconst(1)}; op = {0: sconst(1)}; qp = {0: sconst(1)}
    for e in range(1, 2*dd+2):
        tp[e] = smul(tp[e-1], t); op[e] = smul(op[e-1], omt)
    for e in range(1, dd+1): qp[e] = smul(qp[e-1], q2)
    S = [GQ(0)]*N
    for j in range(1, dd+1):
        inner = [GQ(0)]*N
        for s_ in range(0, j+1):
            c = comb(j, s_)*w(s_, mm)*w(j-s_, dd-mm)
            if c == 0: continue
            inner = sadd(inner, sscale(smul(tp[2*(j-s_)], op[2*s_]), c))
        S = sadd(S, sscale(smul(qp[dd-j], inner), binom_half(j)))
    vp = spow(v, 2*dd-1)
    vpi = [I_ + GQ(0)] ; vi = sadd(v, sconst(I_))        # v + i = 2i + s
    vi = [GQ(0, 2), GQ(1)] + [GQ(0)]*(N-2)
    den = sadd(sadd(spow(v, 4), sscale(spow(v, 2), -6)), sconst(1))
    A = smul(smul(smul(vp, S), spow(vi, 2-2*dd)), sinv(den))
    return [A[n]/A[0] for n in range(N)]
# compile closed forms
a_rf = {dl: RF(T3.arho(dl)) for dl in range(1, 7)}
bad = []; cells = {}
for dd in range(5, 13):
    for mm in range(0, dd+1):
        got = rho_ratios(dd, mm)
        for dl in range(1, 7):
            cells[dl] = cells.get(dl, 0)+1
            if not (got[dl] == a_rf[dl](dd, mm)): bad.append((dd, mm, dl))
print(f'first-principles Taylor series of A(v) at v=i vs closed forms a_delta (delta=1..6): cells {dict(sorted(cells.items()))}, bad {len(bad)}')
print('first bad:', bad[:8], round(time.time()-t0, 1), 's')
