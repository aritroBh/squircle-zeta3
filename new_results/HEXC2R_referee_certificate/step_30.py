
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import reflib as R
lower = {}
for n in (1,2):
    rows, counts = R.make_relations(n, lower)
    lower[n] = rows
print("no-cmzv6 weight2:", len(rows), counts)
