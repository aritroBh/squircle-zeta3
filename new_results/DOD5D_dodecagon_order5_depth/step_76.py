import sys
sys.path.insert(0, ".")
from qs3_symbol import *
L0 = ZERO; L1 = ONE; P0 = q3(4,2); P1 = q3(-3,-2); P2 = q3(F(-1,2),F(-1,2)); P3 = q3(F(3,2),F(1,2))
s1 = g_symbol(L0, (L0, P0), L1)
ref = scale(li_symbol(2, q3div(ONE, P0)), MONE)
print("S G(0,p0;1) = -S Li_2(1/p0):", add(s1, scale(ref, MONE)) == {})
s3 = g_symbol(L0, (L0,)*4+(P0,), L1)
ref3 = scale(li_symbol(5, q3div(ONE, P0)), MONE)
print("S G(0^4,p0;1) = -S Li_5(1/p0):", add(s3, scale(ref3, MONE)) == {})
s2a = shuffle(g_symbol(L0, (P0,), L1), g_symbol(L0, (P2,), L1))
s2b = add(g_symbol(L0, (P0, P2), L1), g_symbol(L0, (P2, P0), L1))
print("shuffle (p0)x(p2) = (p0,p2)+(p2,p0):", add(s2a, scale(s2b, MONE)) == {})
# scaling invariance on a sample: S G_0(w;1) vs S G_0(lambda w; lambda) -- skip, tested below numerically instead
# check the itab/gtab progress
print("--- itab ---"); print(open("itab_build.log").read())
print("--- gtab ---"); print(open("gtab_build.log").read())
