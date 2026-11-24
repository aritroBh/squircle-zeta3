import sys, json, itertools
sys.path.insert(0, '<project>/research/08_theory/GOLD5_golden_corner_w5')
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 40
r5 = mp.sqrt(5)

adata = json.load(open('atom_cobrackets_golden.json'))
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
OBS = {k: tuple((Fr(o[0]), Fr(o[1]) if len(o) > 1 else Fr(0)) for o in adata[k]['obs']) for k in keys}

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

def nfmul(a, b):  # NF mult in Q(sqrt5)
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
    print(f"==== {tag}: {len(sols)} pair solutions, min height {sols[0][0]} ====")
    seen = set()
    for h, k1, k2, s in sols:
        if h > 3 and tag == 'gol1': break
        if h > 9 and tag == 'gol3': break
        d1 = (s[0], s[1]); d2 = (s[2], s[3])
        # verify obs(N5) == obs target
        o = [nfmul(d1, OBS[k1][0]), nfmul(d2, OBS[k2][0])]
        oQ = (o[0][0] + o[1][0], o[0][1] + o[1][1])
        o = [nfmul(d1, OBS[k1][1]), nfmul(d2, OBS[k2][1])]
        oT = (o[0][0] + o[1][0], o[0][1] + o[1][1])
        mtch = (oQ == obs_target[tag][0] and oT == obs_target[tag][1])
        def ds(d): return f"({d[0]}{'+' if d[1]>=0 else ''}{d[1]}r5)"
        def short(k):
            x, y = adata[k]['xy']; return f"Li32({x},{y})"
        key = (k1, k2)
        if key in seen: continue
        seen.add(key)
        print(f"  h={h} obs_match={mtch}: {ds(d1)}*{short(k1)} + {ds(d2)}*{short(k2)}")
