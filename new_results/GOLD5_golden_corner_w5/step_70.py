import sys, json, itertools
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G

adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
y1 = [k for k in keys if adata[k]['xy'][1] == '(Fraction(1, 1), Fraction(0, 1))']
print(f"{len(y1)} atoms with y=1")
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
def cols_of(k):
    v1, v2 = V[k]
    return [(v1, Fr(0), v2, Fr(0)), (Fr(0), v1, Fr(0), v2)]
def height(sol):
    return max(max(abs(x.numerator), x.denominator) for x in sol if x)
T1 = (Fr(-72), Fr(0), Fr(72), Fr(-54))
sols = []
for k1, k2 in itertools.combinations(y1, 2):
    s = solve4(cols_of(k1) + cols_of(k2), T1)
    if s is None: continue
    sols.append((height(s), k1, k2, s))
sols.sort(key=lambda t: t[0])
print(f"gol1 y=1-atom pair solutions: {len(sols)}, best heights: {[s[0] for s in sols[:8]]}")
for h, k1, k2, s in sols[:4]:
    def ds(d): return f"({d[0]}{'+' if d[1]>=0 else ''}{d[1]}r5)"
    print(f"  h={h}: {ds((s[0],s[1]))}*{adata[k1]['xy'][0]},1 + {ds((s[2],s[3]))}*{adata[k2]['xy'][0]},1")
