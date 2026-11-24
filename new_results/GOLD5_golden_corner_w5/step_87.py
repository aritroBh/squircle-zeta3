import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
atomcl = {k: tuple((Fr(q[0]), Fr(q[1])) for q in v)
          for k, v in json.load(open(GOLD + 'atom_classes_golden.json')).items()}
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfsub(a, b):
    return (a[0]-b[0], a[1]-b[1])
SQ5v = (Fr(0), Fr(1))
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
def cols_z(z):
    cl = li5cl[z]
    return [q12(cl), q12(tuple(nfmul(SQ5v, q) for q in cl))]
cl_om = ((Fr(24), Fr(-24)), (Fr(72), Fr(0)), (Fr(144), Fr(0)),
         (Fr(72), Fr(-54)), (Fr(159), Fr(-111)), (Fr(129), Fr(21)))
def solve4(cols, T):
    A = [[cols[c][r] for c in range(4)] + [T[r]] for r in range(4)]
    for col in range(4):
        p = next((r for r in range(col, 4) if A[r][col]), None)
        if p is None: return None
        A[col], A[p] = A[p], A[col]
        q = A[col][col]; A[col] = [x / q for x in A[col]]
        for r in range(4):
            if r != col and A[r][col]:
                A[r] = [a - b * A[r][col] for a, b in zip(A[r], A[col])]
    return [A[i][4] for i in range(4)]
def solve_lin(cols, T):
    n, m = len(T), len(cols)
    A = [[cols[c][r] for c in range(m)] + [T[r]] for r in range(n)]
    R = 0; piv = []
    for c in range(m):
        p = next((r for r in range(R, n) if A[r][c]), None)
        if p is None: continue
        A[R], A[p] = A[p], A[R]
        q = A[R][c]; A[R] = [x / q for x in A[R]]
        for r in range(n):
            if r != R and A[r][c]:
                A[r] = [a - b * A[R][c] for a, b in zip(A[r], A[R])]
        piv.append(c); R += 1
        if R == n: break
    for r in range(R, n):
        if A[r][m]: return None
    part = [Fr(0)] * m
    for i, c in enumerate(piv): part[c] = A[i][m]
    return part
def cols_of(k):
    v1, v2 = V[k]
    return [(v1, Fr(0), v2, Fr(0)), (Fr(0), v1, Fr(0), v2)]
def height(sol):
    return max(max(abs(x.numerator), x.denominator) for x in sol if x)
B4 = [((Fr(3,4),Fr(-1,4)), (Fr(1,4),Fr(-1,4)), (Fr(1,4),Fr(1,4)), (Fr(-1,1),Fr(-1,1))),
      ((Fr(3,4),Fr(-1,4)), (Fr(1,4),Fr(-1,4)), (Fr(1,4),Fr(1,4)), (Fr(-1,1),Fr(1,1))),
      ((Fr(3,4),Fr(-1,4)), (Fr(1,4),Fr(-1,4)), (Fr(-1,1),Fr(1,1)), (Fr(-1,1),Fr(-1,1)))]
B4cols = []
for combo in B4:
    cols = []
    for z in combo: cols += cols_z(z)
    B4cols.append((combo, cols))
T1 = (Fr(-72), Fr(0), Fr(72), Fr(-54))
sols = []
for k1, k2 in itertools.combinations(keys, 2):
    s = solve4(cols_of(k1) + cols_of(k2), T1)
    if s is None: continue
    d1 = (s[0], s[1]); d2 = (s[2], s[3])
    Tcl = tuple(nfsub(q1, nfadd2) for q1, nfadd2 in zip(
        cl_om, tuple((a[0]+b[0], a[1]+b[1]) for a, b in zip(
            tuple(nfmul(d1, q) for q in atomcl[k1]), tuple(nfmul(d2, q) for q in atomcl[k2])))))
    T12 = q12(Tcl)
    best_h = None
    for combo, cols in B4cols:
        sp = solve_lin(cols, T12)
        if sp is None: continue
        h = height(sp)
        if best_h is None or h < best_h[0]:
            best_h = (h, combo, sp)
    hn5 = max(max(abs(x.numerator), x.denominator) for x in s if x)
    sols.append((best_h[0] if best_h else 10**18, hn5, k1, k2, (d1, d2), best_h))
sols.sort(key=lambda t: (t[0], t[1]))
print(f"{len(sols)} gol1 N5-pairs; best by P5 height:")
for h5, hn5, k1, k2, (d1, d2), bh in sols[:6]:
    print(f"  P5h={h5} N5h={hn5}: ({G.nf_str(d1)})*Li32{tuple(adata[k1]['xy'])} + ({G.nf_str(d2)})*Li32{tuple(adata[k2]['xy'])}")
    if bh:
        h, combo, sp = bh
        for i, z in enumerate(combo):
            dd = (sp[2*i], sp[2*i+1])
            if any(dd): print(f"       ({G.nf_str(dd)})*Li5({z})")
