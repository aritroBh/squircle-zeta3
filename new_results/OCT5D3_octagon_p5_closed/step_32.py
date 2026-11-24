import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
from oct5d3_basis import *
d = pickle.load(open("symfit_m4f.pkl","rb"))
x0, K, labels = d['x0'], d['K'], d['vis_labels']
nz = [(labels[i], x0[i]) for i in range(len(x0)) if x0[i] != ZERO]
print("RREF particular solution: %d nonzero terms" % len(nz))
maxnum = 1
for lab, v in nz:
    maxnum = max(maxnum, abs(v.p.numerator), abs(v.q.numerator))
    print("  %-26s %s" % (lab, v))
print("max |numerator|:", maxnum)
# kernel vector heights
print("\nkernel vectors (first 5) with heights:")
from fractions import Fraction
for v in K[:5]:
    nzs = [(labels[i], v[i]) for i in range(len(v)) if v[i]!=ZERO]
    h = max([1]+[abs(z.p.numerator) for _,z in nzs]+[abs(z.q.numerator) for _,z in nzs])
    print("  height", h, ":", [(l, str(z)) for l,z in nzs][:6])
