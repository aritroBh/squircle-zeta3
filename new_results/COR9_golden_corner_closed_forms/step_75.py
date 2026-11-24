from fractions import Fraction as Fr
vecs = [(0, 2, 0, 0), (-2, 0, 1, -2), (-1, 2, 0, -1), (-1, 5, 0, -1)]
# cor7_pt's letter_span_rank elimination
M = [list(map(Fr, v)) for v in vecs]
R = 0
for c in range(len(M[0])):
    p_ = next((r for r in range(R, len(M)) if M[r][c]), None)
    if p_ is None: continue
    M[R], M[p_] = M[p_], M[R]
    q = M[R][c]; M[R] = [x/q for x in M[R]]
    for r in range(len(M)):
        if r != R and M[r][c]:
            M[r] = [a - b*M[R][c] for a, b in zip(M[r], M[R])]
    R += 1
print("cor7_pt-style rank:", R)
# independent check
import sympy as sp
print("sympy rank:", sp.Matrix([list(v) for v in vecs]).rank())
