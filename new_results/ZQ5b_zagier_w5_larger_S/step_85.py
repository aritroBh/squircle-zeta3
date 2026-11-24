import os, json
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
import sys
sys.path.insert(0,'.')
from w5blib import QD
J = json.load(open('out/d13/candidates.json'))
cands=[QD(Ff(c['x'][0]),Ff(c['x'][1]),13) for c in J['cands']]
S=set(cands)
# x is duplication-redundant if x = y^2 with y in S and -y in S (then Li5(x)=16(Li5(y)+Li5(-y)))
redundant=set()
for i,x in enumerate(cands):
    # find square roots: y with y^2==x.  Check candidates y
    for y in S:
        if y*y==x and -y in S:
            redundant.add(i); break
print(f"total {len(cands)}, duplication-redundant (x=y^2, +-y in pool): {len(redundant)}")
# also inversion-related? keep simple. How many remain
print(f"remaining after removing redundant: {len(cands)-len(redundant)}")
# what about x such that x=y^2 but only one of +-y in pool? still expressible? Li5(x)=16(Li5(y)+Li5(-y)) needs both.
# count pairs
npairs=0
for i,x in enumerate(cands):
    for y in S:
        if y*y==x:
            npairs+=1
print("total square relations x=y^2 in pool:", npairs)
