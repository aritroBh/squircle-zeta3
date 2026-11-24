import sys; sys.path.insert(0, '<project>/research/08_theory/COR9_golden_corner_closed_forms')
import importlib, cor9_rank; importlib.reload(cor9_rank)
from cor9_rank import letter_vecs, rank_int, kernel_basis, hensel_mu
from fractions import Fraction as Fr
import mpmath as mp
D, a, b, e = 2, 1, 1, 3
vecs, primes = letter_vecs(D, a, b, e)
pids = sorted({pid for v in vecs for pid in v})
M = [[v.get(pid, 0) for pid in pids] for v in vecs]
print("pids:", pids)
for row in M: print(row)
print("rank:", rank_int(M), " kernel:", kernel_basis(M))
