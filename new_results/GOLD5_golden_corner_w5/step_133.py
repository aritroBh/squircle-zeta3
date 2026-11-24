import sys, itertools, json, pickle
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G

def is_lyndon(w):
    return all(w < w[i:] + w[:i] for i in range(1, len(w)))
LW = {n: [w for w in itertools.product((0, 1), repeat=n) if is_lyndon(w)] for n in range(1, 6)}
ALL = sorted([w for n in range(1, 6) for w in LW[n]])
def monoms(n):
    out = []
    for r in range(1, n + 1):
        for combo in itertools.combinations_with_replacement(ALL, r):
            if sum(len(w) for w in combo) == n:
                out.append(tuple(sorted(combo, reverse=True)))
    return out
M5 = monoms(5)
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
MEXP = {m: expand_mon(m) for m in M5}
lead = {tuple(x for L in m for x in L): m for m in M5}
def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfsub(a, b): return (a[0]-b[0], a[1]-b[1])
def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
Z2 = (Fr(0), Fr(0))
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
    assert not rem
    return out
W3M = [m for m in M5 if any(len(L) == 3 for L in m)]
def w3vec(SE):
    lex = lyndon_expand(SE)
    return tuple(pad(lex.get(m, Z2)) for m in W3M)
print("W3M monomials:", [(''.join(map(str,L)) for L in m) for m in W3M])
P = G.golden_point('gol3')
w3_om = w3vec(G.ft_expand_sym_nf(P.omega5(), 2))
print("Omega_5 gol3 w3-content:", [G.nf_str(q) for q in w3_om])
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
w3_atoms = {}
for k in keys:
    x, y = eval(adata[k]['xy'][0]), eval(adata[k]['xy'][1])
    w3_atoms[k] = w3vec(G.ft_expand_sym_nf(G.li32_sym_nf(P, x, y), 2))
print("atom w3-contents done")
li5cl = {eval(k): None for k in json.load(open(GOLD + 'li5_classes_golden.json'))}
import mpmath as mp
mp.mp.dps = 30
def fval(z):
    return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))
nzsmall = [z for z in li5cl if abs(fval(z)) <= 1]
w3_li5 = {}
for z in nzsmall:
    w3_li5[z] = w3vec(G.ft_expand_sym_nf(G.li5_sym_nf(P, z), 2))
nzw3 = [z for z, v in w3_li5.items() if any(any(q) for q in v)]
print("Li5 args with nonzero w3-content:", [(str(z), fval(z)) for z in nzw3])
pickle.dump({'w3_om': w3_om, 'w3_atoms': w3_atoms, 'w3_li5': w3_li5, 'W3M': W3M},
            open(GOLD + 'w3_content_gol3.pkl', 'wb'))
