"""check_rho_def.py -- independent re-derivation of the initial data rho_{pq}(delta) (level k=d) from the
integrand definition: S_m(t) polynomial from WTS genintegrand.S_poly (a different implementation of S than
TOWN1's trunc_engine), Common(s) = -(sqrt2/2) 8^d (i+s)^(2d-1) (2i+s)^(2-2d)/G(i+s), t(v)=(v^2+2v-1)/(4v),
G=v^4-6v^2+1.  rho(delta) = [s^delta] Common(s) [S_p(t) + [p>q] S_q(1-t)].
Compared (exactly) with (i) TOWN1 pickle level k=d, (ii) the TOW3 closed forms a_5, a_6.
"""
import sys, pickle, time
from fractions import Fraction as Fr
sys.path.insert(0, '<project>/research/08_theory/WTS_weight_k_structure/theory')
import genintegrand as GI   # read-only import (-B: no bytecode)

D = 6
class G:
    __slots__ = ('r', 'i')
    def __init__(s, r=0, i=0): s.r = Fr(r); s.i = Fr(i)
    def __add__(s, o): return G(s.r+o.r, s.i+o.i)
    def __sub__(s, o): return G(s.r-o.r, s.i-o.i)
    def __mul__(s, o): return G(s.r*o.r - s.i*o.i, s.r*o.i + s.i*o.r)
    def inv(s):
        n = s.r*s.r + s.i*s.i; return G(s.r/n, -s.i/n)
    def iszero(s): return s.r == 0 and s.i == 0
def smul(a, b):
    out = [G() for _ in range(D+1)]
    for x, ax in enumerate(a):
        if ax.iszero(): continue
        for y, by in enumerate(b):
            if x+y > D: break
            out[x+y] = out[x+y] + ax*by
    return out
def sinv(a):
    a0i = a[0].inv(); out = [G() for _ in range(D+1)]
    out[0] = a0i
    for n in range(1, D+1):
        acc = G()
        for k in range(1, n+1): acc = acc + a[k]*out[n-k]
        out[n] = (G()-acc)*a0i
    return out
def sconst(c): return [G(c)] + [G() for _ in range(D)]
def spow(a, n):
    out = sconst(1)
    for _ in range(n): out = smul(out, a)
    return out
def sadd(a, b): return [x+y for x, y in zip(a, b)]
def sscale(a, c): return [x*G(c) for x in a]

ONE = G(1)
def ser(*cs):
    return [c if isinstance(c, G) else G(c) for c in cs] + [G() for _ in range(D+1-len(cs))]
v = ser(G(0, 1), 1)                       # v = i + s
v2 = smul(v, v)
numer = sadd(sadd(v2, sscale(v, 2)), sconst(-1))   # v^2 + 2v - 1
t = smul(numer, sinv(sscale(v, 4)))       # t(v)
omt = sadd(sconst(1), sscale(t, -1))      # 1 - t
Gv = sadd(sadd(smul(v2, v2), sscale(v2, -6)), sconst(1))
Ginv = sinv(Gv)

def horner(coeffs, x):
    out = sconst(0)
    for c in reversed(coeffs):
        out = sadd(smul(out, x), sconst(c))
    return out

def rho_indep(d, p, q):
    Sp = GI.S_poly(d, p)
    ser_ = horner(Sp, t)
    if p > q:
        Sq = GI.S_poly(d, q)
        ser_ = sadd(ser_, horner(Sq, omt))
    com = smul(smul(spow(v, 2*d-1), sinv(spow(ser(G(0, 2), 1), 2*d-2))), Ginv)
    X = smul(com, ser_)
    pref = Fr(-(8**d), 2)         # Common has factor -(sqrt2/2) 8^d = sqrt2 * pref
    return [(Fr(0), pref*x.r, Fr(0), pref*x.i) for x in X]    # (p0,p1,p2,p3) of sqrt2*pref*X

if __name__ == '__main__':
    t0 = time.time()
    H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl', 'rb'))
    import check_data as CD   # reuse Q class + closed-form parser
    a5, a6 = CD.load_md_closed_forms()
    ev = {5: CD.make_eval(a5, (CD.d_, CD.m_)), 6: CD.make_eval(a6, (CD.d_, CD.m_))}
    tot = 0; bad = 0; badCF = 0; totCF = 0
    for d in range(5, 13):
        for p in range(d, (d-1)//2, -1):
            q = d-p
            r = rho_indep(d, p, q)
            for dl in range(0, min(D, 2*d-4)+1):
                tot += 1
                dat = tuple(H[(p, q)][dl][d])
                if tuple(r[dl]) != dat:
                    bad += 1
                    if bad <= 5: print('BAD', d, p, q, dl, r[dl], dat)
            r0 = CD.Q(*r[0])
            for dl in (5, 6):
                rr = CD.Q(*r[dl]) * r0.inv()
                er, ei = ev[dl](d, p)
                totCF += 1
                if not (rr.in_Qi() and rr.a == er and rr.c == ei):
                    badCF += 1; print('CF-BAD', d, p, q, dl)
    print(f'independent rho(delta<=6) from S_poly-definition vs TOWN1 pickle (level k=d), d=5..12: {tot} coefficients, bad = {bad}')
    print(f'independent rho(5)/rho(0), rho(6)/rho(0) vs proof_depth6.md closed forms, d=5..12: {totCF} cells, bad = {badCF}')
    print('time', round(time.time()-t0, 1))
