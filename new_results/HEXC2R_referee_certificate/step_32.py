
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import reflib as R, cmzv6
from fractions import Fraction as Fr

extra = {(((1, 0), 1), ((2, 0), 1), ((3, 0), 1), ((4, 0), 1), ((5, 0), 1)),
         (((3, 0), 1),), (((-1, 3), 1), ((3, 3), 1))}
for d in (2,3,6):
    mine = [R.primitive(r) for r in R.distribution(2, d)]
    theirs = [cmzv6.primitive(r) for r in cmzv6.distribution(2, d)]
    print("d =", d, "mine", len([r for r in mine if r]), "theirs", len([r for r in theirs if r]))
    sm, st = set(mine), set(theirs)
    if sm != st:
        for r in sm - st: print("  M-only:", r)
        for r in st - sm: print("  T-only:", r)
