import sys, json, itertools
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr

adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
def nf_of(o):
    if not o: return (Fr(0), Fr(0))
    return (Fr(o[0]), Fr(o[1]) if len(o) > 1 else Fr(0))
OBS = {k: (nf_of(adata[k]['obs'][0]), nf_of(adata[k]['obs'][1])) for k in keys}

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

def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])

obs_target = {'gol1': ((Fr(72), Fr(0)), (Fr(72), Fr(-54))),
              'gol3': ((Fr(0), Fr(0)), (Fr(972), Fr(432)))}
targets = {'gol1': (Fr(-72), Fr(0), Fr(72), Fr(-54)),
           'gol3': (Fr(0), Fr(0), Fr(972), Fr(432))}
for tag, T in targets.items():
    sols = []
    for k1, k2 in itertools.combinations(keys, 2):
        s = solve4(cols_of(k1) + cols_of(k2), T)
        if s is None: continue
        sols.append((height(s), k1, k2, s))
    sols.sort(key=lambda t: t[0])
    hmin = sols[0][0]
    print(f"==== {tag}: {len(sols)} pair solutions, min height {hmin} ====")
    for h, k1, k2, s in sols:
        if h > hmin: break
        d1 = (s[0], s[1]); d2 = (s[2], s[3])
        oQ = tuple(sum(t) for t in zip(nfmul(d1, OBS[k1][0]), nfmul(d2, OBS[k2][0])))
        oT = tuple(sum(t) for t in zip(nfmul(d1, OBS[k1][1]), nfmul(d2, OBS[k2][1])))
        mtch = (oQ == obs_target[tag][0] and oT == obs_target[tag][1])
        def ds(d): return f"({d[0]}{'+' if d[1]>=0 else ''}{d[1]}r5)"
        def short(k):
            x, y = adata[k]['xy']; return f"Li32({x},{y})"
        print(f"  h={h} obs_match={mtch}: {ds(d1)}*{short(k1)} + {ds(d2)}*{short(k2)}")
