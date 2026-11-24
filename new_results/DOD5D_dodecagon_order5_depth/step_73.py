import sys, random
sys.path.insert(0, ".")
from fractions import Fraction as F
from qs3_symbol import *
import obstruction32 as OB

# --- control: Li_5 at Q(sqrt3) S-units pairs to zero (r=3) ---
pts = [q3(2,-1), q3(F(-1,2),F(1,2)), q3(0,F(1,3)), q3(7,-4), q3(1,F(-1,2)),
       q3(-1,F(2,3)), q3(F(1,2),F(-1,6)), q3(F(1,3)), q3(1), q3(F(1,2)),
       q3(-1,1), q3(3,2), q3(5,2)]
ok = True
for z in pts:
    try:
        s = li_symbol(5, z)
    except Exception as e:
        print("logvec FAIL at", z, ":", e); ok = False; continue
    table, _ = OB.obstruction(s, 3)
    if not OB.is_zero(table):
        print("NONZERO obstruction for Li_5 at", z); ok = False
print("Li_5 at", len(pts), "S-units pairs to 0:", ok)

# --- control: random shuffle products pair to zero ---
random.seed(11)
okB = True
for _ in range(300):
    k = random.choice([1, 2])
    w1 = tuple(random.choice([0,1,2]) for _ in range(k))
    w2 = tuple(random.choice([0,1,2]) for _ in range(5-k))
    s1 = {w1: q3(random.randint(-4,4), random.randint(-3,3))}
    s2 = {w2: q3(random.randint(-4,4), random.randint(-3,3))}
    table, _ = OB.obstruction(shuffle(s1, s2), 3)
    if not OB.is_zero(table): okB = False; break
print("300 random shuffle products pair to 0:", okB)

# --- G-symbol identities on test words (12-gon letters) ---
L0 = ZERO; L1 = ONE; P0 = q3(4,2); P1 = q3(-3,-2)
s1 = g_symbol(L0, (L0, P0), L1)
ref = scale(li_symbol(2, q3div(ONE, P0)), MONE)
print("S G(0,p0;1) = -S Li_2(1/p0):", add(s1, scale(ref, MONE)) == {})
s3 = g_symbol(L0, (L0,)*4+(P0,), L1)
ref3 = scale(li_symbol(5, q3div(ONE, P0)), MONE)
print("S G(0^4,p0;1) = -S Li_5(1/p0):", add(s3, scale(ref3, MONE)) == {})
s2a = shuffle(g_symbol(L0, (P0,), L1), g_symbol(L0, (P1,), L1))
s2b = add(g_symbol(L0, (P0, P1), L1), g_symbol(L0, (P1, P0), L1))
print("shuffle (p0)x(p1) = (p0,p1)+(p1,p0):", add(s2a, scale(s2b, MONE)) == {})
