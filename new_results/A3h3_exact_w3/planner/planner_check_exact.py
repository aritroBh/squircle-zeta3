"""Planner re-check of A3h3 with independent code: exact Q(sqrt2)-ranks of E2 and E2+R_i via the regular
representation [[A,2B],[B,A]] and flint fmpz_mat.rank. Also a toy test of the block formula."""
import pickle
from fractions import Fraction as F
from math import lcm
from flint import fmpz_mat
T = 'research/08_theory/'
w = pickle.load(open(T + 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb'))
k = pickle.load(open(T + 'A3f3_level8_compute/k3_newHrows.pickle', 'rb'))['krows']
R = pickle.load(open(T + 'A3f4_glanois_basis/five_relations.pickle', 'rb'))
cols = list(w['cols']); idx = {c: i for i, c in enumerate(cols)}; n = len(cols)

def block_rank(rows):  # rows: list of dict col -> (a, b) meaning a + b*sqrt2
    M = []
    for r in rows:
        den = lcm(*[F(x).denominator for ab in r.values() for x in ab]) if r else 1
        top = [0] * (2 * n); bot = [0] * (2 * n)
        for c, (a, b) in r.items():
            a = int(F(a) * den); b = int(F(b) * den)
            top[c] += a; top[n + c] += 2 * b      # row (A | 2B)
            bot[c] += b; bot[n + c] += a          # row (B | A)
        M.append(top); M.append(bot)
    return fmpz_mat(M).rank() // 2

# toy: rows (1, sqrt2), (sqrt2, 2) -> rank 1 over Q(sqrt2)
n_save = n; n = 2
assert block_rank([{0: (1, 0), 1: (0, 1)}, {0: (0, 1), 1: (2, 0)}]) == 1
assert block_rank([{0: (1, 0), 1: (0, 1)}, {0: (0, 1), 1: (3, 0)}]) == 2
n = n_save
E1 = [{idx[wd]: (c, 0) for wd, c in r} for r in w['rows']]
E2 = E1 + [{idx[wd]: (x, y) for wd, (x, y) in r.items()} for r in k]
# R_i as vectors on the 32 words (7,b,c), b outer 4..7, c inner 0..7
print('R structure:', type(R).__name__, (list(R.keys())[:5] if isinstance(R, dict) else len(R)))
vecs = R['relations'] if isinstance(R, dict) and 'relations' in R else (R if isinstance(R, list) else list(R.values())[0])
order = [(7, b, c) for b in range(4, 8) for c in range(8)]
Rrows = [{idx[wd]: (int(v), 0) for wd, v in zip(order, vec) if v} for vec in vecs]
r2 = block_rank(E2); print('exact rank E2 =', r2)
print('exact rank E2 + R_i =', [block_rank(E2 + [Ri]) for Ri in Rrows])
print('exact rank E2 + all R =', block_rank(E2 + Rrows))
