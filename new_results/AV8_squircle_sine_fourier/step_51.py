from mpmath import mp, mpf, pi, log, nstr
mp.dps = 60
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2)
# corrected: (denom, a, b, c) with b_n'' = (a + b ln2 + c pi^2)/denom
rels2 = {3:(128,12,-16,3), 5:(3456,508,48,-45), 7:(13824,-572,-96,63),
         9:(64000,2532,160,-225), 11:(1296000,-25592,-2160,2475), 13:(88905600,1598344,75600,-143325)}
from fractions import Fraction
from mpmath import rational
print("k=2 relation check:")
for n,(d,a,b,c) in rels2.items():
    v = mpf(res[f"{n},2"])
    pred = (a + b*ln2 + c*pi**2)/d
    g = Fraction(rational.gcd(rational.gcd(a,b),rational.gcd(b,c)))
    print(n, nstr(v-pred,4), "  simplified:", (a+g-a)//g, (b)//g, (c)//g, "/", d//g)
