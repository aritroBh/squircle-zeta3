
import json
from fractions import Fraction as F
import dec5d_lib as D
from qsqrt2 import QS2

kit2 = D.make_kit(QS2, D.logvec2, ('log2','log(1+sqrt2)'))
obs2, L2_2, L3_2 = D.make_obstruction(kit2)
ZERO2 = QS2(0)
symj = json.loads(open("<project>/research/08_theory/M4c_w5_closed_form/symbol.json").read())['symbol']
Sm4 = {tuple(map(int,w)): QS2(F(p), F(q)) for w,(p,q) in symj.items()}
obs_m4 = obs2(Sm4)
print("M4c W5 obstruction with DEC5D general code:")
print("  vs [X,[X,Y]]x[X,Y]:", obs_m4[0][0])
print("  vs [Y,[X,Y]]x[X,Y]:", obs_m4[1][0])
print("expected: 84 and -38*sqrt2:", QS2(84), QS2(0,-38))
print("PASS" if obs_m4[0][0]==QS2(84) and obs_m4[1][0]==QS2(0,-38) and len(L2_2)==1 and len(L3_2)==2 else "FAIL")
# also Li_5 and shuffle controls over Q(sqrt2)
zs2 = [QS2(-1,1), -QS2(-1,1), -QS2(-1,1)**2, D.SQ2/2, QS2(F(1,2)), 1-D.SQ2/2, 2-D.SQ2]
okA = all(all(v == ZERO2 for row in obs2(kit2['li_symbol'](5, z)) for v in row) for z in zs2)
print("control Li_5 at 7 Q(sqrt2) S-units:", "PASS" if okA else "FAIL")
import random
random.seed(11); okB = True
for _ in range(300):
    k = random.choice([1,2])
    w1 = tuple(random.choice([0,1]) for _ in range(k))
    w2 = tuple(random.choice([0,1]) for _ in range(5-k))
    s1 = {w1: QS2(random.randint(-4,4), random.randint(-3,3))}
    s2 = {w2: QS2(random.randint(-4,4), random.randint(-3,3))}
    if not all(v == ZERO2 for row in obs2(kit2['shuffle'](s1,s2)) for v in row): okB = False; break
print("control 300 shuffle products over Q(sqrt2):", "PASS" if okB else "FAIL")
