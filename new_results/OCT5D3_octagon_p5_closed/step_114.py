import os, sys, pickle, math
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
from fractions import Fraction as F
rel1500 = {k: tuple(F(x) for x in v) for k, v in pickle.load(open("rel1500_frac.pkl","rb")).items()}
def lcm(a,b): return a*b//math.gcd(a,b)
L = 1
maxnum = 1
for k,(p,q) in rel1500.items():
    L = lcm(L, p.denominator); L = lcm(L, q.denominator)
    maxnum = max(maxnum, abs(p.numerator), abs(q.numerator))
print("closure relation: %d terms; lcm of denominators = %d; max |numerator| = %d" % (len(rel1500), L, maxnum))
print("lcm factorization:", end=" ")
n = L
facs = {}
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]:
    while n % p == 0:
        facs[p] = facs.get(p,0)+1; n //= p
print(facs, "rest:", n)
# the LLL integer relation height (target coeff v0): recover v0 = lcm? The fit reported maxc=2.2e9.
# terms
print("\nLchi5 coefficient:", rel1500['Lchi5'])
print("z5 coefficient:", rel1500['z5'])
