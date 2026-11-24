
import mpmath as mp, sys, time
mp.mp.dps = 110
sys.path.insert(0, "research/08_theory/D17_zeta5_polygons")
import poly_coeffs5 as p5
print("dps:", mp.mp.dps)
m = 3
t0=time.time()
pc = p5.regular_P_coeffs5(m, 5)   # [d^k]P = 4m int_0^{A/2} D_k
from math import factorial
P5 = 120*pc[5]
print("P5(m=3) =", mp.nstr(P5, 60), " (%.0fs)"%(time.time()-t0))
# closed-form predictions
h = mp.pi/(2*m); Sh=mp.sin(h); Ch=mp.cos(h)
C5 = 2880*m*Sh**4*Ch*(65*Sh**6-63*Sh**4+17*Sh**2-1)
C23 = -960*m*Ch*Sh**4*(45*Sh**6-27*Sh**4+6*Sh**2-1)
print("C_3^(5) =", mp.nstr(C5, 30), "  2835sqrt3/32 =", mp.nstr(2835*mp.sqrt(3)/32,30))
print("C_3^(2,3) =", mp.nstr(C23,30), "  1395sqrt3/32 =", mp.nstr(1395*mp.sqrt(3)/32,30))
rem = P5 - C5*mp.zeta(5) - C23*mp.zeta(2)*mp.zeta(3)
print("rem =", mp.nstr(rem, 60))
import json
json.dump({"P5":mp.nstr(P5,100),"C5":mp.nstr(C5,100),"C23":mp.nstr(C23,100),"rem":mp.nstr(rem,100)}, open("m3_total.json","w"))
