import sys, itertools, json, time
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
import mpmath as mp
mp.mp.dps = 30
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
atomcl = {k: tuple((Fr(q[0]), Fr(q[1])) for q in v)
          for k, v in json.load(open(GOLD + 'atom_classes_golden.json')).items()}
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
def fval(z):
    return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))
nzsmall = [z for z, cl in li5cl.items() if abs(fval(z)) <= 1 and any(any(q) for q in cl)]
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
SQ5v = (Fr(0), Fr(1))
def cols_z(z):
    cl = li5cl[z]
    return [q12(cl), q12(tuple(nfmul(SQ5v, q) for q in cl))]
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
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[R])]
        piv.append(c); R += 1
        if R == n: break
    for r in range(R, n):
        if A[r][m]: return None
    part = [Fr(0)] * m
    for i, c in enumerate(piv): part[c] = A[i][m]
    return part
def solve4(cols, T):
    A = [[cols[c][r] for c in range(4)] + [T[r]] for r in range(4)]
    for col in range(4):
        p = next((r for r in range(col, 4) if A[r][col]), None)
        if p is None: return None
        A[col], A[p] = A[p], A[col]
        q = A[col][col]; A[col] = [x / q for x in A[col]]
        for r in range(4):
            if r != col and A[r][col]:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return [A[i][4] for i in range(4)]
def cols_of(k):
    v1, v2 = V[k]
    return [(v1, Fr(0), v2, Fr(0)), (Fr(0), v1, Fr(0), v2)]
def height(sol):
    return max(max(abs(x.numerator), x.denominator) for x in sol if x)
bases = []
for combo in itertools.combinations(nzsmall, 4):
    cols = []
    for z in combo: cols += cols_z(z)
    bases.append((combo, cols))
cl_om = {'gol1': ((Fr(24), Fr(-24)), (Fr(72), Fr(0)), (Fr(144), Fr(0)),
                  (Fr(72), Fr(-54)), (Fr(159), Fr(-111)), (Fr(129), Fr(21))),
         'gol3': ((Fr(-6), Fr(0)), (Fr(0), Fr(0)), (Fr(-21), Fr(3)), (Fr(972), Fr(432)),
                  (Fr(-84), Fr(18)), (Fr(-10563), Fr(-4449)))}
T_del = {'gol1': (Fr(-72), Fr(0), Fr(72), Fr(-54)),
         'gol3': (Fr(0), Fr(0), Fr(972), Fr(432))}
t0 = time.time()
for tag in ('gol1', 'gol3'):
    sols = []
    for k1, k2 in itertools.combinations(keys, 2):
        s = solve4(cols_of(k1) + cols_of(k2), T_del[tag])
        if s is None: continue
        d1 = (s[0], s[1]); d2 = (s[2], s[3])
        Tcl = tuple((a[0]-b1[0]-b2[0], a[1]-b1[1]-b2[1]) for a, b1, b2 in zip(
            cl_om[tag], tuple(nfmul(d1, q) for q in atomcl[k1]), tuple(nfmul(d2, q) for q in atomcl[k2])))
        T12 = q12(Tcl)
        bh = None
        for combo, cols in bases:
            sp = solve_lin(cols, T12)
            if sp is None: continue
            h = height(sp)
            if bh is None or h < bh[0]: bh = (h, combo, sp)
        hn5 = max(max(abs(x.numerator), x.denominator) for x in s if x)
        sols.append((bh[0] if bh else 10**18, hn5, k1, k2, (d1, d2), bh))
    sols.sort(key=lambda t: (t[0] + t[1], t[0]))
    print(f"==== {tag} ({time.time()-t0:.0f}s): best by P5h+N5h:")
    for h5, hn5, k1, k2, (d1, d2), bh in sols[:5]:
        print(f"  P5h={h5} N5h={hn5}: ({G.nf_str(d1)})*Li32{tuple(adata[k1]['xy'])} + ({G.nf_str(d2)})*Li32{tuple(adata[k2]['xy'])}")
        if bh:
            h, combo, sp = bh
            for i, z in enumerate(combo):
                dd = (sp[2*i], sp[2*i+1])
                if any(dd): print(f"       ({G.nf_str(dd)})*Li5({z})")
