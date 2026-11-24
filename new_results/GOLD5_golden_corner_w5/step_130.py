import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G
from cor6_sym import sshuffle

def is_lyndon(w):
    return all(w < w[i:] + w[:i] for i in range(1, len(w)))
LW = {n: [w for w in itertools.product((0, 1), repeat=n) if is_lyndon(w)] for n in range(1, 6)}
def monoms(n):
    out = []
    def rec(rem, start, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        for k in range(min(rem, start), 0, -1):
            for L in LW[k]:
                if cur and L > cur[-1]: continue
                rec(rem - k, k, cur + [L])
    rec(n, 5, [])
    return out
M5 = monoms(5)
print(f"{len(M5)} monomials (expect 32)")
def expand_mon(m):
    S = {(): Fr(1)}
    for L in m:
        S = sshuffle(S, {L: Fr(1)})
    return S
# triangularity check: leading word of monomial = concatenation in DECREASING order
lead = {}
ok = True
for m in M5:
    cat = tuple(x for L in m for x in L)
    E = expand_mon(m)
    mx = max(E)
    lead[cat] = m
    if mx != cat: ok = False; print("nontriangular:", m, ''.join(map(str, mx)))
print("triangular (leading = decreasing concat):", ok)
# invert: express raw symbol in monomial basis
def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfsub(a, b): return (a[0]-b[0], a[1]-b[1])
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
Z2 = (Fr(0), Fr(0))
MEXP = {m: expand_mon(m) for m in M5}
def lyndon_expand(SE):
    rem = {w: pad(cf) for w, cf in SE.items() if any(pad(cf))}
    out = {}
    for lw in sorted(lead, reverse=True):
        if lw not in rem: continue
        cf = rem.pop(lw)
        m = lead[lw]
        out[m] = cf
        for w, c in MEXP[m].items():
            if w == lw: continue
            if w in rem:
                rem[w] = nfsub(rem[w], nfmul(cf, (c, Fr(0))))
                if not any(rem[w]): del rem[w]
    assert not rem, f"leftover: {list(rem)[:3]}"
    return out
MRHO  = (Fr(1,2), Fr(-1,2)); MPH2  = (Fr(1,4), Fr(-1,4))
PH2   = (Fr(1,4), Fr(1,4)); TP2   = (Fr(3), Fr(-1)); TP2O2 = (Fr(3,4), Fr(-1,4))
RHO3 = (Fr(-2), Fr(1)); RHO6 = (Fr(9), Fr(-4)); M1 = (Fr(-1), Fr(0)); ONE = (Fr(1), Fr(0))
N5 = {'gol1': [((Fr(-8, 3), Fr(0)), MRHO, MPH2), ((Fr(8, 27), Fr(2, 3)), M1, RHO3)],
      'gol3': [((Fr(-9, 8), Fr(-1, 2)), RHO6, ONE)]}
P5 = {'gol1': [((Fr(25, 12), Fr(-37, 12)), TP2O2), ((Fr(-17), Fr(9)), MPH2),
               ((Fr(16, 3), Fr(-4, 3)), PH2), ((Fr(-19, 12), Fr(-3, 4)), TP2)],
      'gol3': [((Fr(11, 6), Fr(-1, 2)), MPH2), ((Fr(17, 2), Fr(-5, 2)), PH2),
               ((Fr(20, 9), Fr(-16, 27)), RHO3), ((Fr(-26, 3), Fr(2)), TP2)]}
for tag in ('gol1', 'gol3'):
    P = G.golden_point(tag)
    S = {}
    for w, cf in P.omega5().items():
        S[w] = nfadd(S.get(w, Z2), pad(cf))
    for d, x, y in N5[tag]:
        for w, cf in G.li32_sym_nf(P, x, y).items():
            S[w] = nfsub(S.get(w, Z2), nfmul(d, pad(cf if isinstance(cf, tuple) else (cf,))))
    for c, z in P5[tag]:
        for w, cf in G.li5_sym_nf(P, z).items():
            S[w] = nfsub(S.get(w, Z2), nfmul(c, pad(cf)))
    SE = G.ft_expand_sym_nf(S, 2)
    lex = lyndon_expand(SE)
    w3 = {m: cf for m, cf in lex.items() if any(len(L) == 3 for L in m)}
    print(f"==== {tag}: {len(lex)} monomials; with weight-3 Lyndon factor: {len(w3)}")
    for m, cf in sorted(w3.items()):
        print(f"   {[''.join(map(str,L)) for L in m]}: {G.nf_str(cf)}")
