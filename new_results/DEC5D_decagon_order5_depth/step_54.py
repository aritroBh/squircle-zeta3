
import dec5d_lib as D
from fractions import Fraction as F
QS5 = D.QS5; PHI = D.PHI
# integrand constants: log(2 phi^2), log(2 phi)
print("logvec5(2 phi^2) =", D.logvec5(2*PHI**2), " logvec5(2 phi) =", D.logvec5(2*PHI))

kit5 = D.make_kit(QS5, D.logvec5, ('log2','logsqrt5','logphi'))
obs3, L2_3, L3_3 = D.make_obstruction(kit5)
print("d=3 obstruction: %d Lie_3 x %d Lie_2 = %d pairings" % (len(L3_3), len(L2_3), len(L3_3)*len(L2_3)))
ZERO = QS5(0)
def obs3_zero(s):
    return all(v == ZERO for row in obs3(s) for v in row)

# control (a): Li_5 at Q(sqrt5) S-units, including sqrt5-involving ones
SQ5 = D.SQ5
zs = [QS5(F(1,2)), QS5(2), -QS5(1), 1/PHI, 1/PHI**2, PHI, PHI**2, -PHI, -1/PHI, SQ5, 1/SQ5, SQ5/5, 2/PHI, (3-SQ5)/2]
ok = True
for z in zs:
    z = QS5.c(z)
    # check 1-z is S-unit (logvec must not raise)
    s = kit5['li_symbol'](5, z)
    good = obs3_zero(s)
    ok &= good
    print(f"  Li_5({z}): pairings all zero: {good}")
print("control Li_5: ", "PASS" if ok else "FAIL")

# control (b): random shuffle products over the appearing slots
import random
random.seed(17)
okB = True
for _ in range(300):
    k = random.choice([1,2,3])
    slots = [0,2] if random.random()<0.7 else [0,1,2]
    w1 = tuple(random.choice(slots) for _ in range(k))
    w2 = tuple(random.choice(slots) for _ in range(5-k))
    s1 = {w1: QS5(random.randint(-4,4), random.randint(-3,3))}
    s2 = {w2: QS5(random.randint(-4,4), random.randint(-3,3))}
    if not obs3_zero(kit5['shuffle'](s1,s2)): okB = False; break
print("control 300 random shuffle products pair to 0:", "PASS" if okB else "FAIL")
