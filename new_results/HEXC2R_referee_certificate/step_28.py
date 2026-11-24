
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import reflib as R
from fractions import Fr if False else __import__("fractions").Fraction as Fr
# fresh: build relations without any prior cache
lower = {}
counts_fresh = {}
for n in (1,2,3):
    rows, counts = R.make_relations(n, lower)
    lower[n] = rows
    counts_fresh[n] = (len(rows), counts)
print("fresh:", counts_fresh)
