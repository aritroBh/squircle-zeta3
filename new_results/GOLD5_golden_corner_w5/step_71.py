import sys, itertools
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 40
r5 = mp.sqrt(5)

P = G.golden_point('gol1')
F = P.F
ONE = F.one()
def fval(x):
    return mp.mpf(x[0].numerator)/x[0].denominator + mp.mpf(x[1].numerator)/x[1].denominator*r5
def vec_or_none(el):
    try: return P.vec(el)
    except ValueError: return None
sols = [ONE, F.neg(ONE)]
for tup in itertools.product(range(-9, 10), repeat=2):
    if tup == (0,0): continue
    v = ONE
    for e_, g in zip(tup, G.BASIS):
        if e_: v = F.mul(v, F.pow(g, e_))
    for sgn in (1, -1):
        x = v if sgn == 1 else F.neg(v)
        if F.eq(x, ONE): continue
        om = F.sub(ONE, x)
        if all(t == 0 for t in om): continue
        if vec_or_none(om) is not None and not any(F.eq(x, y) for y in sols):
            sols.append(x)
zs = [s for s in sols if not F.eq(s, ONE)]
zc = {}
for z in zs:
    S = G.li3_sym_nf(P, z)
    SE = G.ft_expand_sym_nf(S, 2)
    qc, _ = G.quotient_coords_nf(SE)
    zc[str(z)] = qc
print("Li3 quotient coords (word (1,0,0), word (1,1,0)) as NF:")
for k, v in sorted(zc.items(), key=lambda kv: float(abs(fval(eval(kv[0]))))):
    print(f"  z={k:26s}: ({G.nf_str(v[0])}), ({G.nf_str(v[1])})")

# solve Delta = sum c_z Li3(z), c_z in F, over Q^4
def q4(v):  # NF^2 -> Q^4
    return (v[0][0], v[0][1], v[1][0], v[1][1])
def solve4(cols, T):
    M = [[cols[c][r] for c in range(4)] for r in range(4)]
    A = [row[:] + [T[i]] for i, row in enumerate(M)]
    for col in range(4):
        p = next((r for r in range(col, 4) if A[r][col]), None)
        if p is None: return None
        A[col], A[p] = A[p], A[col]
        q = A[col][col]; A[col] = [x / q for x in A[col]]
        for r in range(4):
            if r != col and A[r][col]:
                A[r] = [a - b * A[r][col] for a, b in zip(A[r], A[col])]
    return [A[i][4] for i in range(4)]
SQ5v = (Fr(0), Fr(1))
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def cols_z(k):
    v = zc[k]
    return [q4(v), q4((nfmul(SQ5v, v[0]), nfmul(SQ5v, v[1])))]
targets = {'gol1': (Fr(-72), Fr(0), Fr(72), Fr(-54)),
           'gol3': (Fr(0), Fr(0), Fr(972), Fr(432))}
zkeys = list(zc.keys())
def height(sol):
    return max(max(abs(x.numerator), x.denominator) for x in sol if x)
for tag, T in targets.items():
    best = []
    for k1, k2 in itertools.combinations(zkeys, 2):
        s = solve4(cols_z(k1) + cols_z(k2), T)
        if s is None: continue
        best.append((height(s), k1, k2, s))
    best.sort(key=lambda t: t[0])
    print(f"==== {tag}: Li3-part solve, {len(best)} pairs, best:")
    for h, k1, k2, s in best[:4]:
        d1 = (s[0], s[1]); d2 = (s[2], s[3])
        print(f"  h={h}: ({G.nf_str(d1)})*Li3({k1}) + ({G.nf_str(d2)})*Li3({k2})")
