import sys, itertools, json
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from fractions import Fraction as Fr
from cor7_pt import QuadField, QuadPoint
from cor6_sym import W4Quotient, lneg

Q2 = W4Quotient(2)
def li4_coords(P, x):
    vx, vo = P.vec(x), P.vec(P.F.sub(P.F.one(), x))
    return Q2.coords({(lneg(vo), vx, vx, vx): Fr(1)})

QF = QuadField(5)
P = QuadPoint(QF, "phi", (Fr(-1,2), Fr(1,2)), primes=[2])
args = {"1/phi": (Fr(-1,2),Fr(1,2)), "phi": (Fr(1,2),Fr(1,2)), "-phi": (Fr(-1,2),Fr(-1,2)),
        "1/phi^2": (Fr(3,2),Fr(-1,2)), "phi^2": (Fr(3,2),Fr(1,2)), "-1/phi": (Fr(1,2),Fr(-1,2)),
        "2": (Fr(2),Fr(0)), "1/2": (Fr(1,2),Fr(0)), "-1": (Fr(-1),Fr(0)),
        "phi/2": (Fr(1,4),Fr(1,4)), "2/phi": (Fr(-1),Fr(1)), "-phi^2": (Fr(-3,2),Fr(-1,2)),
        "1/(2phi^2)": (Fr(1,4),Fr(-1,4))}
cls = {nm: li4_coords(P, x) for nm, x in args.items()}
target = [(Fr(111,2), Fr(0), Fr(-12)), (Fr(-51,2), Fr(0), Fr(12))]
def solve3(trip, v):
    M = [[cls[n][j] for n in trip] for j in range(3)]
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
found = None
for trip in itertools.combinations(cls, 3):
    s0 = solve3(trip, list(target[0]))
    if s0 is None: continue
    s1 = solve3(trip, list(target[1]))
    if s1 is None: continue
    score = sum(abs(x.numerator)+x.denominator for x in s0+s1)
    if found is None or score < found[0]:
        found = (score, trip, s0, s1)
score, trip, s0, s1 = found
print("golden corner: class(Omega(1/phi)) components = Li4 classes:")
print("  theta^0:", " + ".join(f"({c})Li4({n})" for c, n in zip(s0, trip) if c))
print("  theta^1:", " + ".join(f"({c})Li4({n})" for c, n in zip(s1, trip) if c))
json.dump({"trip": trip, "s0": [str(x) for x in s0], "s1": [str(x) for x in s1]},
          open("golden_li4_combo.json", "w"))
