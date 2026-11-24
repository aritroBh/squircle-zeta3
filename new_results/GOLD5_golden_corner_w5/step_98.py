import sys, itertools, json
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
def height(sol):
    return max(max(abs(x.numerator), x.denominator) for x in sol if x)
keyof = {tuple(v['xy']): k for k, v in adata.items()}
cl_om = {'gol1': ((Fr(24), Fr(-24)), (Fr(72), Fr(0)), (Fr(144), Fr(0)),
                  (Fr(72), Fr(-54)), (Fr(159), Fr(-111)), (Fr(129), Fr(21))),
         'gol3': ((Fr(-6), Fr(0)), (Fr(0), Fr(0)), (Fr(-21), Fr(3)), (Fr(972), Fr(432)),
                  (Fr(-84), Fr(18)), (Fr(-10563), Fr(-4449)))}
N5c = {
 'gol1': [((Fr(-1,4), Fr(0)), '(Fraction(1, 2), Fraction(-1, 4))', '(Fraction(1, 1), Fraction(0, 1))'),
          ((Fr(4,9), Fr(2,3)), '(Fraction(-1, 1), Fraction(0, 1))', '(Fraction(-2, 1), Fraction(1, 1))')],
 'gol3': [((Fr(-9,8), Fr(-1,2)), '(Fraction(9, 1), Fraction(-4, 1))', '(Fraction(1, 1), Fraction(0, 1))')],
}
for tag in ('gol1', 'gol3'):
    Tcl = list(cl_om[tag])
    for d, xs, ys in N5c[tag]:
        k = keyof[(xs, ys)]
        Tcl = [(a[0]-b[0], a[1]-b[1]) for a, b in zip(Tcl, [nfmul(d, q) for q in atomcl[k]])]
    T12 = q12(tuple(Tcl))
    best = []
    for r in (4, 5, 6):
        for combo in itertools.combinations(nzsmall, r):
            cols = []
            for z in combo: cols += cols_z(z)
            sp = solve_lin(cols, T12)
            if sp is None: continue
            # VERIFY
            acc = [Fr(0)]*12
            for j, c in enumerate(sp):
                if c:
                    for rr in range(12): acc[rr] += c*cols[j][rr]
            assert all(x == 0 for x in acc)
            best.append((height(sp), combo, sp))
    best.sort(key=lambda t: t[0])
    print(f"==== {tag}: target {[G.nf_str(q) for q in Tcl]}")
    print(f"   verified solutions: {len(best)}, best heights {[b[0] for b in best[:8]]}")
    for h, combo, sp in best[:2]:
        print(f"  h={h}:")
        for i, z in enumerate(combo):
            dd = (sp[2*i], sp[2*i+1])
            if any(dd): print(f"    ({G.nf_str(dd)})*Li5({z})")
