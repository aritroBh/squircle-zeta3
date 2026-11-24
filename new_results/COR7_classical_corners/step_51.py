import sys
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L
from fractions import Fraction as Fr

mp.mp.dps = 1050
W3 = mp.mpf(open("c13_W3_hpl_1150.txt").read())
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
basis3_names = ["zeta3", "Li3(1/3)", "Li3(-1/2)", "pi2*log2", "pi2*log3",
                "Li2(1/3)*log2", "Li2(1/3)*log3", "log2^3", "log2^2*log3", "log2*log3^2", "log3^3"]
basis3 = [mp.zeta(3), Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2), mp.pi**2*l2, mp.pi**2*l3,
          Li(2, mp.mpf(1)/3)*l2, Li(2, mp.mpf(1)/3)*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]

cands = {}
for dps in (600, 900):
    rows = L.find_relation_rows([W3] + basis3, dps)
    cands[dps] = L.normalize(rows[0])
    mp.mp.dps = 1050
    resid = sum(c*v for c, v in zip(cands[dps], [W3] + basis3))
    print(f"dps {dps}: height {max(abs(c) for c in cands[dps])}, resid@1050 {mp.nstr(abs(resid),3)}")
print("identical at 600/900:", cands[600] == cands[900])
a = cands[900]; a0 = a[0]
print("\nW3 =")
for i in range(1, len(a)):
    if a[i]:
        print(f"  {str(Fr(-a[i], a0)):>8} * {basis3_names[i-1]}")
