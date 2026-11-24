
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import cmzv6
import reflib as R
lower = {}
for n in (1,2):
    rows, counts = R.make_relations(n, lower)
    lower[n] = rows
print("with-cmzv6 weight2:", len(rows), counts)
lowert = {}
for n in (1,2):
    rowst, countst = cmzv6.make_relations(n, lowert)
    lowert[n] = rowst
print("cmzv6 weight2:", len(rowst), countst)
print("mine - theirs:", set(rows)-set(rowst))
