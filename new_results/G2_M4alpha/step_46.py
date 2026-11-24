import time
from fractions import Fraction as Fr
import flint, random
# benchmark fmpq rref at 521 x 485 with realistic-size rationals
random.seed(1)
rows, cols = 521, 486
M = flint.fmpq_mat(rows, cols)
for i in range(rows):
    for j in range(cols):
        M[i, j] = flint.fmpq(random.randint(-50, 50), random.randint(1, 200))
t0 = time.time()
R, rank = M.rref()
print("521x486 rref:", time.time() - t0, "rank", rank)
