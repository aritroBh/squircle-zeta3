import sys, itertools, json, pickle
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
Fraction = Fr
import gold5_lib as G
from cor6_sym import rref
import mpmath as mp
mp.mp.dps = 30

def solve_packed(pklpath, cl_om, co_om):
    D = pickle.load(open(pklpath, 'rb'))
    w3_om, w3_atoms, w3_li5 = D['w3_om'], D['w3_atoms'], D['w3_li5']
    adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
    atomcl = {k: tuple((Fr(q[0]), Fr(q[1])) for q in v)
              for k, v in json.load(open(GOLD + 'atom_classes_golden.json')).items()}
    li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
             for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
    keys = list(w3_atoms.keys())
    def fval(z):
        return float(mp.mpf(z[0].numerator)/z[0].denominator + mp.mpf(z[1].numerator)/z[1].denominator*mp.sqrt(5))
    nzsmall = [z for z in li5cl if abs(fval(z)) <= 1]
    V = {k: (Fr(adata[k]['cvec'][0][0]), Fr(adata[k]['cvec'][1][0])) for k in keys}
    def pad(a): return a if len(a) == 2 else (a[0], Fr(0))
    def nfmul(a, b):
        return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
    SQ5v = (Fr(0), Fr(1))
    def pack_atom(k):
        nfv = [(V[k][0], Fr(0)), (V[k][1], Fr(0))]
        nfv += [pad(q) for q in atomcl[k]]
        nfv += [pad(q) for q in w3_atoms[k]]
        return nfv
    def pack_li5(z):
        nfv = [(Fr(0), Fr(0)), (Fr(0), Fr(0))]
        nfv += [pad(q) for q in li5cl[z]]
        nfv += [pad(q) for q in w3_li5[z]]
        return nfv
    Tom = [co_om] + [pad(q) for q in cl_om] + [pad(q) for q in w3_om]
    def qvec(nfv):
        out = []
        for a in nfv: out += [a[0], a[1]]
        return tuple(out)
    Tq = qvec(Tom)
    cols = []; names = []
    for k in keys:
        nfv = pack_atom(k)
        cols.append(qvec(nfv)); names.append(('atom', k, 0))
        cols.append(qvec(tuple(nfmul(SQ5v, a) for a in nfv))); names.append(('atom', k, 1))
    for z in nzsmall:
        nfv = pack_li5(z)
        cols.append(qvec(nfv)); names.append(('li5', z, 0))
        cols.append(qvec(tuple(nfmul(SQ5v, a) for a in nfv))); names.append(('li5', z, 1))
    return cols, names, Tq

def solve_lin_full(cols, T):
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
        if A[r][m]: return None, None
    part = [Fr(0)] * m
    for i, c in enumerate(piv): part[c] = A[i][m]
    free = [c for c in range(m) if c not in piv]
    ker = []
    for f in free:
        v = [Fr(0)] * m; v[f] = Fr(1)
        for i, c in enumerate(piv): v[c] = -A[i][f]
        ker.append(v)
    return part, ker
def height(sol):
    return max((max(abs(x.numerator), x.denominator) for x in sol if x), default=0)

for tag, pkl, cl_om, co_om in [
    ('gol1', GOLD + 'w3_content_gol1.pkl',
     ((Fr(24), Fr(-24)), (Fr(72), Fr(0)), (Fr(144), Fr(0)), (Fr(72), Fr(-54)), (Fr(159), Fr(-111)), (Fr(129), Fr(21))),
     (Fr(-72), Fr(0)), ),
    ('gol3', GOLD + 'w3_content_gol3.pkl',
     ((Fr(-6), Fr(0)), (Fr(0), Fr(0)), (Fr(-21), Fr(3)), (Fr(972), Fr(432)), (Fr(-84), Fr(18)), (Fr(-10563), Fr(-4449))),
     (Fr(0), Fr(0)), )]:
    pass
# need cobracket pair: gol1 = ((-72,0),(72,-54)); gol3 = ((0,0),(972,432))
for tag, pkl, cl_om, co_om in [
    ('gol1', GOLD + 'w3_content_gol1.pkl',
     ((Fr(24), Fr(-24)), (Fr(72), Fr(0)), (Fr(144), Fr(0)), (Fr(72), Fr(-54)), (Fr(159), Fr(-111)), (Fr(129), Fr(21))),
     ((Fr(-72), Fr(0)), (Fr(72), Fr(-54)))),
    ('gol3', GOLD + 'w3_content_gol3.pkl',
     ((Fr(-6), Fr(0)), (Fr(0), Fr(0)), (Fr(-21), Fr(3)), (Fr(972), Fr(432)), (Fr(-84), Fr(18)), (Fr(-10563), Fr(-4449))),
     ((Fr(0), Fr(0)), (Fr(972), Fr(432))))]:
    cols, names, Tq = solve_packed(pkl, cl_om, co_om)
    part, ker = solve_lin_full(cols, Tq)
    if part is None:
        print(f"{tag}: NOT solvable"); continue
    # verify
    acc = [Fr(0)]*len(Tq)
    for j, c in enumerate(part):
        if c:
            for rr in range(len(Tq)): acc[rr] += c*cols[j][rr]
    ok = all(a-b == 0 for a, b in zip(acc, Tq))
    print(f"==== {tag}: pivot-particular solution (verify {ok}), height {height(part)}, kernel dim {len(ker)}")
    terms = [(names[j], part[j]) for j in range(len(part)) if part[j]]
    print(f"  {len(terms)} nonzero terms:")
    for nm, c in terms:
        kind, ident, ab = nm
        coef = (c, Fr(0)) if ab == 0 else (Fr(0), c)
        if kind == 'atom':
            print(f"    ({G.nf_str(coef)})*Li32{tuple(json.load(open(GOLD + 'atom_cobrackets_golden.json'))[ident]['xy'])}")
        else:
            print(f"    ({G.nf_str(coef)})*Li5({ident})")
