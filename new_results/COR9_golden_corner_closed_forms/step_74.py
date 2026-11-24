import sys
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
from fractions import Fraction as Fr
from cor9_pt import CFQuadField, QuadPoint
from cor9_rank import true_rank, letter_vecs, rank_int, kernel_basis
import mpmath as mp
D, a, b, e = 2, 9, -4, 7
# by hand: true rank
r, M, pids, rel = true_rank(D, a, b, e)
print("true_rank:", r)
for row, pid in zip(M, ["2", "c", "1-c", "1+c"]): print("  ", pid, row, pids)
# cor7_pt
QF = CFQuadField(D)
P = QuadPoint(QF, "x", (Fr(a, e), Fr(b, e)), primes=[2, 7])
print("gensE:", [(str(g[0]), g[1]) for g in P.gensE])
rr, labs, vecs = P.letter_span_rank()
print("cor7_pt rank:", rr)
for nm, v in vecs: print("  ", nm, v)
