"""W4RX: independent numerical evaluation of all level-8 values G(a1,...,an;1), n<=4.

Own notation (used in every W4RX script): a word is a tuple of letters; letter m in 0..7 is the
point zeta^m (zeta=exp(2 pi i/8)), letter 8 is the point 0.
    G(a1,...,an; y) = int_0^y dt/(t-a1) G(a2,...,an; t),  G(;y)=1.
A word is admissible (convergent at both ends) iff a1 != zeta^0 = 1 and an != 0.

Method (no code of the notes is used): Hoelder convolution at p=2,
    G(a1..an;1) = sum_k (-1)^k G(1-ak,...,1-a1; 1/2) G(a(k+1),...,an; 1/2),
valid for a1 != 1, an != 0, and the power series in y of G(b1..bm;y) (bm != 0), whose
coefficients obey, for a prepended letter c:
    c = 0 :  g[n] = g_w[n]/n
    c != 0:  h[n] = (h[n-1] + g_w[n-1])/c,   g[n] = -h[n]/n.
All letters c that occur satisfy |c| >= |1-zeta| = 0.765, so the series at y=1/2 converge like
0.6533^n.  Arithmetic: python-flint acb.  Run twice (N terms and N2 terms) to measure truncation.
Usage: python evalG.py N PREC_BITS out.json
"""
import sys, json, time
from itertools import product
from flint import acb, arb, ctx

N = int(sys.argv[1]); ctx.prec = int(sys.argv[2]); OUT = sys.argv[3]
ZERO = 8
def zeta(m):
    m %= 8
    exact = {0: acb(1), 2: acb(0, 1), 4: acb(-1), 6: acb(0, -1)}
    return exact[m] if m in exact else acb(0, arb.pi() * m / 4).exp()
Z = [zeta(m) for m in range(8)]
Y = acb(1) / 2
YPOW = [acb(1)]
for n in range(1, N + 1): YPOW.append(YPOW[-1] * Y)
INV = [None] + [acb(1) / n for n in range(1, N + 1)]

def prepend(c, g):
    """Series of G(c,w;y) from the series g of G(w;y) (g has length N+1)."""
    new = [acb(0)] * (N + 1)
    if c == 0:
        assert g[0] == 0
        for n in range(1, N + 1): new[n] = g[n] * INV[n]
        return new
    ic = 1 / c; h = acb(0)
    for n in range(1, N + 1):
        h = (h + g[n - 1]) * ic
        new[n] = -h * INV[n]
    return new

def value(g):
    s = acb(0)
    for n in range(N + 1): s += g[n] * YPOW[n]
    return s


t0 = time.time()
# Right pieces: R[w] = G(w; 1/2), w over letters 0..8, last letter != 8.  Built by prepending.
R = {(): acb(1)}
def dfs_right(word, g):
    if len(word) == 4: return
    for c in range(9):
        new = (c,) + word
        gn = prepend(acb(0) if c == ZERO else Z[c], g)
        R[new] = value(gn)
        dfs_right(new, gn)
g0 = [acb(1)] + [acb(0)] * N
for c in range(8):          # last letter must be nonzero
    gn = prepend(Z[c], g0); R[(c,)] = value(gn); dfs_right((c,), gn)
print(f'right pieces {len(R)} {time.time()-t0:.1f}s', flush=True)
# Left pieces: L[(a1..ak)] = G(1-ak,...,1-a1; 1/2), a1 != zeta^0.  Prepending 1-a(k+1).
L = {(): acb(1)}
def one_minus(c): return acb(1) if c == ZERO else 1 - Z[c]
def dfs_left(prefix, g):
    if len(prefix) == 4: return
    for c in range(9):
        new = prefix + (c,)
        gn = prepend(one_minus(c) if not (c == 0) else acb(0), g)  # 1 - zeta^0 = 0
        L[new] = value(gn)
        dfs_left(new, gn)
for c in range(1, 9):       # a1 != 1
    gn = prepend(one_minus(c), g0); L[(c,)] = value(gn); dfs_left((c,), gn)
print(f'left pieces {len(L)} {time.time()-t0:.1f}s', flush=True)

def G1(w):
    s = acb(0)
    for k in range(len(w) + 1):
        term = L[w[:k]] * R[w[k:]]
        s += -term if k % 2 else term
    return s

out = {}
for n in range(1, 5):
    for w in product(range(9), repeat=n):
        if w[0] == 0 or w[-1] == ZERO: continue
        v = G1(w)
        out[' '.join(map(str, w))] = [v.real.mid().str(170, radius=False), v.imag.mid().str(170, radius=False),
                                      v.real.rad().str(3, radius=False), v.imag.rad().str(3, radius=False)]
json.dump({'N': N, 'prec': ctx.prec, 'values': out}, open(OUT, 'w'))
print(f'done words={len(out)} {time.time()-t0:.1f}s', flush=True)
