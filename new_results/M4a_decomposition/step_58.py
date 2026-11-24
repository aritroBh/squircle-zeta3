import os
BASE = "<project>"
p = os.path.join(BASE, "research/08_theory/A3_W4_proof/gpl_expression.tsv")
lines = open(p).read().split("\n")
print("\n".join(lines[:8]))
print("...")
print("\n".join(lines[-4:]))
print("total lines:", len(lines))
