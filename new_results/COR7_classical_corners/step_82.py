import sys, json
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 60

# ---- angles of the classical corners ----
corners = [("1/2 (hexagon)", mp.mpf(1)/2), ("1/3", mp.mpf(1)/3),
           ("sqrt2-1", mp.sqrt(2)-1), ("(sqrt2-1)^2", (mp.sqrt(2)-1)**2),
           ("1/sqrt2 (square)", 1/mp.sqrt(2)), ("1/phi (golden)", (mp.sqrt(5)-1)/2),
           ("1/phi^3 = sqrt5-2", mp.sqrt(5)-2)]
print("=== classical corners and angles ===")
for nm, c in corners:
    a = 2*mp.acos(c)
    print(f"c0 = {nm:18s} c={mp.nstr(c,15)}  alpha = {mp.nstr(a,30)} rad = {mp.nstr(a*180/mp.pi,15)} deg")

# ---- explicit Li4 combinations for the golden corner components ----
from cor7_pt import QuadField, QuadPoint
from cor6_sym import W4Quotient, lneg
Q2 = W4Quotient(2)
def li4_coords(P, x):
    vx, vo = P.vec(x), P.vec(P.F.sub(P.F.one(), x))
    return Q2.coords({(lneg(vo), vx, vx, vx): Fr(1)})
QF = QuadField(5)
P = QuadPoint(QF, "phi", (Fr(-1,2), Fr(1,2)), primes=[2])
F = P.F
args = {"1/phi": (Fr(-1,2),Fr(1,2)), "phi": (Fr(1,2),Fr(1,2)), "-phi": (Fr(-1,2),Fr(-1,2)),
        "1/phi^2": (Fr(3,2),Fr(-1,2)), "phi^2": (Fr(3,2),Fr(1,2)), "-1/phi": (Fr(1,2),Fr(-1,2)),
        "2": (Fr(2),Fr(0)), "1/2": (Fr(1,2),Fr(0)), "-1": (Fr(-1),Fr(0)),
        "phi/2": (Fr(1,4),Fr(1,4)), "2/phi": (Fr(1),Fr(1)), "-phi^2": (Fr(-3,2),Fr(-1,2))}
cls = {nm: li4_coords(P, x) for nm, x in args.items()}
import itertools
target = [(Fr(111,2), Fr(0), Fr(-12)), (Fr(-51,2), Fr(0), Fr(12))]
# solve with 3 args
names = list(cls)
def solve3(trip, v):
    M = [[cls[n][j] for n in trip] for j in range(3)]
    # solve M c = v
    A = [row[:] + [v[i]] for i, row in enumerate(M)]
    for col in range(3):
        p = next((r for r in range(col,3) if A[r][col]), None)
        if p is None: return None
        A[col], A[p] = A[p], A[col]
        q = A[col][col]; A[col] = [x/q for x in A[col]]
        for r in range(3):
            if r != col and A[r][col]:
                A[r] = [x - y*A[r][col] for x, y in zip(A[r], A[col])]
    return [A[i][3] for i in range(3)]
for trip in itertools.combinations(names, 3):
    s0 = solve3(trip, list(target[0]))
    if s0 is None: continue
    s1 = solve3(trip, list(target[1]))
    if s1 is None: continue
    if all(abs(x.numerator) < 2000 and x.denominator < 200 for x in s0+s1):
        print("\ngolden corner: class(Omega) = class of")
        print("  theta^0:", " + ".join(f"({c})Li4({n})" for c, n in zip(s0, trip) if c))
        print("  theta^1:", " + ".join(f"({c})Li4({n})" for c, n in zip(s1, trip) if c))
        break

# ---- W2 at sqrt2-1: exact formula + candidate reductions ----
mp.mp.dps = 120
c = mp.sqrt(2)-1
W2v = mp.mpf(open("b21_W2_hpl_1150.txt").read()[:130])
W2f = -mp.polylog(2, c) - mp.polylog(2, -c) + mp.pi**2/12
print("\nW2(sqrt2-1): |stored - (-Li2(b)-Li2(-b)+pi^2/12)| =", mp.nstr(abs(W2v-W2f), 3))
# candidate: Li2(b) + Li2(-b) = 1/2 Li2(b^2)? duplication: Li2(x)+Li2(-x) = 1/2 Li2(x^2)
lhs = mp.polylog(2, c) + mp.polylog(2, -c)
rhs = mp.polylog(2, c*c)/2
print("duplication Li2(b)+Li2(-b) - 1/2 Li2(b^2) =", mp.nstr(abs(lhs-rhs), 3))
# so W2 = -1/2 Li2((sqrt2-1)^2) + pi^2/12 exactly (classical duplication)
