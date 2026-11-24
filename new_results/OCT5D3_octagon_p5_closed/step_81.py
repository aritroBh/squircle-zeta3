import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
from qsqrt2 import QS2 as Q
from fractions import Fraction as F
from flint import fmpz_mat
import math

sf = pickle.load(open("symfit_m4f.pkl","rb"))
K, labels = sf['K'], sf['vis_labels']

def lcm(a, b): return a*b//math.gcd(a, b)
# integer embedding: each K vector (52 QS2) -> Z^104 with common denom
ints = []
for v in K:
    D = 1
    for z in v:
        D = lcm(D, z.p.denominator); D = lcm(D, z.q.denominator)
    iv = []
    for z in v:
        iv += [int(z.p*D), int(z.q*D)]
    ints.append((D, iv))
print("denominators:", [D for D,_ in ints])
M = fmpz_mat([iv for _, iv in ints])
Mr = M.lll()
Kred = []
for i in range(Mr.nrows()):
    row = [int(Mr[i,j]) for j in range(Mr.ncols())]
    if any(row): Kred.append(row)
print("reduced kernel vectors:", len(Kred))
for row in Kred:
    h = max(abs(x) for x in row)
    nzs = [(labels[j//2], (row[j] if j%2==0 else 0)) for j in range(104) if row[j]!=0]
    # group by atom
    from collections import defaultdict
    d = defaultdict(list)
    for j in range(104):
        if row[j]: d[labels[j//2]].append(row[j])
    print("  height %-6d terms %d: %s" % (h, len(d), [(k, tuple(v)) for k, v in list(d.items())][:4]))
pickle.dump({'Kred': Kred, 'labels': labels}, open("kernel_red.pkl","wb"))
print("saved kernel_red.pkl")
