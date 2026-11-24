import itertools
def is_lyndon(w):
    return all(w < w[i:] + w[:i] for i in range(1, len(w)))
LW = {n: [w for w in itertools.product((0, 1), repeat=n) if is_lyndon(w)] for n in range(1, 6)}
ALL = sorted([w for n in range(1, 6) for w in LW[n]])
def monoms(n):
    out = []
    for r in range(1, n + 1):
        for combo in itertools.combinations_with_replacement(ALL, r):
            if sum(len(w) for w in combo) == n:
                out.append(combo)
    return out
M5 = monoms(5)
print(f"{len(M5)} monomials (expect 32)")
from fractions import Fraction as Fr
def sshuf(u, v):
    if not u: return [v]
    if not v: return [u]
    return ([(u[0],) + w for w in sshuf(u[1:], v)] + [(v[0],) + w for w in sshuf(u, v[1:])])
def expand_mon(m):
    S = {(): Fr(1)}
    for L in m:
        N = {}
        for u, cu in S.items():
            for w in sshuf(u, L):
                N[w] = N.get(w, Fr(0)) + cu
        S = N
    return S
lead = {}
ok = True
# convention: sort each monomial in DECREASING order for the leading-term property
M5s = []
for m in M5:
    md = tuple(sorted(m, reverse=True))
    M5s.append(md)
    cat = tuple(x for L in md for x in L)
    E = expand_mon(md)
    if max(E) != cat:
        ok = False
    lead[cat] = md
print("triangular:", ok, "unique leads:", len(lead) == len(M5s))
