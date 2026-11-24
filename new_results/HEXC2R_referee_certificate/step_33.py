
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import reflib as R, cmzv6
from fractions import Fraction as Fr
from itertools import product

n, d = 2, 6
N = 6
w = (0, 0)
# my code trace
C = {(j,): Fr(1) for j in range(N // d, N, N // d)}
cp = [{(): Fr(1)}]
for k in range(1, n + 1):
    cp.append(R.vadd((R.mul(cp[-1], C), Fr(1, k))))
print("C =", C)
print("cp[1] =", cp[1])
print("cp[2] =", cp[2])
lift_choices = []
for j in w:
    if j == -1: lift_choices.append(((-1,),))
    else: lift_choices.append(tuple(((j // d + k * (N // d)) % N,) for k in range(d)))
left = {}
for tup in product(*lift_choices):
    w2 = tuple(x[0] for x in tup)
    left[w2] = left.get(w2, Fr(0)) + d ** w.count(-1)
print("left raw:", left)
right = R.reg(w)
for k in range(1, n + 1):
    if w[k - 1] != 0: break
    right = R.vadd((right, 1), (R.mul(cp[k], R.reg(w[k:])), 1))
print("right:", right)
row = R.regularize(R.vadd((left, 1), (right, -1)))
print("row:", row)
print("primitive:", R.primitive(row))
# theirs
rows_t = list(cmzv6.distribution(2, 6))
for r in rows_t:
    print("theirs row:", cmzv6.primitive(r))
