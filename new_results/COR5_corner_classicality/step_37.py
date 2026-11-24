import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import W4Quotient, li_sym, rref
sols = json.load(open("sunit23_solutions.json"))
def vec23(fac):
    s, a, b = fac
    return (Fr(a), Fr(b))
Q = W4Quotient(2)
li4 = [(Fr(xs), Q.coords(li_sym(4, vec23(fx), vec23(fy)))) for xs, fx, fy in sols]
vF = (Fr(5,2), Fr(1,2), Fr(0))
M = [list(v) for _, v in li4] + [list(vF)]
print("input rows:", len(M), "cols:", len(M[0]))
# instrumented rref
rows = [list(r) for r in M]; R = 0; piv = []
for j in range(len(rows[0])):
    p = next((i for i in range(R, len(rows)) if rows[i][j]), None)
    print(f"col {j}: R={R}, pivot row={p}")
    if p is None: continue
    rows[R], rows[p] = rows[p], rows[R]
    q = rows[R][j]; rows[R] = [x/q for x in rows[R]]
    for i in range(len(rows)):
        if i != R and rows[i][j]:
            rows[i] = [a - b*rows[i][j] for a, b in zip(rows[i], rows[R])]
    piv.append(j); R += 1
    if R == len(rows): break
print("final R:", R)
