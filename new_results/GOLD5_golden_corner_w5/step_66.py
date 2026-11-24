import sys, json, itertools
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from fractions import Fraction as Fr

adata = json.load(open('atom_cobrackets_golden.json'))
# rational cvecs
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: tuple(Fr(c) for q in adata[k]['cvec'] for c in q[:1]) for k in keys}  # (v1, v2) rational
for k in keys:
    assert all(Fr(adata[k]['cvec'][i][1]) == 0 for i in range(2)), "cvec not rational"
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}

def solve4(cols, T):
    """solve M x = T, M 4x4 (columns), exact; None if singular."""
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

targets = {'gol1': (Fr(-72), Fr(0), Fr(72), Fr(-54)),
           'gol3': (Fr(0), Fr(0), Fr(972), Fr(432))}
for tag, T in targets.items():
    sols = []
    for k1, k2 in itertools.combinations(keys, 2):
        cols = cols_of(k1) + cols_of(k2)
        s = solve4(cols, T)
        if s is None: continue
        # d1 = a+b sqrt5, d2 = c+d sqrt5
        h = height(s)
        sols.append((h, k1, k2, s))
    sols.sort(key=lambda t: t[0])
    print(f"==== {tag}: target {T} ====")
    for h, k1, k2, s in sols[:6]:
        d1 = f"({s[0]}{'+' if s[1]>=0 else ''}{s[1]}*sqrt5)"
        d2 = f"({s[2]}{'+' if s[3]>=0 else ''}{s[3]}*sqrt5)"
        print(f"  h={h}: N5 = {d1}*{k1} + {d2}*{k2}")
