
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import reflib as R
from fractions import Fraction as Fr
lower = {}
res = {}
for n in (1,2,3):
    rows, counts = R.make_relations(n, lower)
    lower[n] = rows
    res[n] = len(rows)
print("fresh-cache run:", res)

# now pollute cache and repeat in-process with fresh 'lower'
w = (2, 0, 0)
R.reg(w)
lower2 = {}
res2 = {}
for n in (1,2,3):
    rows2, counts2 = R.make_relations(n, lower2)
    lower2[n] = rows2
    res2[n] = len(rows2)
print("after reg((2,0,0)):", res2)
