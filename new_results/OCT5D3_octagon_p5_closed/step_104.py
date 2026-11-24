import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
ak = pickle.load(open("atom_kernel.pkl","rb"))
rels = ak['rels']
print("total candidate rows saved:", len(rels))
import math
for h, v in rels:
    print("h=%.3e" % h if h else 0, end="  ")
print()
# group by height
hs = sorted(h for h, _ in rels)
print("heights:", ["%.1e" % h for h in hs])
