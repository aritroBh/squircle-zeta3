
from fractions import Fraction as F
import engine5, mzv_exact as M

truth = {3:[F(10)], 4:[F(-19)], 5:[F(102),F(-36)], 6:[F(-159),F(68)]}
m0 = m0_coords((0,1,1,1,1), 6)
def peval5(p, a):
    return sum(c*a[0]**e[0]*a[1]**e[1]*a[2]**e[2]*a[3]**e[3]*a[4]**e[4] for e,c in p.items())
forms6, eng6 = engine5.reduced_forms(6)
for W in range(3,7):
    reg = [peval5(forms6[W][b], (0,1,1,1,1)) for b in range(len(truth[W]))]
    m0v = m0[W]
    treg = [truth[W][b]-m0v[b] for b in range(len(truth[W]))]
    print(f"W={W}: truth={truth[W]} m0={m0v} => true regular={treg}  my regular={reg}")
# also print Areg word dicts at W=4,5
for W in [4,5]:
    print(f"--- Areg[{W}] ---")
    for w,p in eng6['Areg'][W].items():
        print("  ", w, p)
