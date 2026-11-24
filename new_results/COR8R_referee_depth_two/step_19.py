from cor8r_num import *
setdps(40)
# Derived formula: with a0 := z, a_{n+1} := 0,
# dG(a1..an;z) = sum_i G(skip a_i; z) dlog((a_{i-1}-a_i)/(a_i-a_{i+1}))
# Check n=3 numerically: G(a,b,cc;z), args chosen generic, vary a,b,cc,z independently.
import random
random.seed(1)
a, b, cc, z = [mp.mpf(v) for v in ('2.7','5.3','9.1','1.0')]
h = mp.mpf('1e-14')
def G3(aa,bb,ccc,zz):
    return G_int([aa,bb,ccc], zz)
# partials
dA = (G3(a+h,b,cc,z)-G3(a-h,b,cc,z))/(2*h)
dB = (G3(a,b+h,cc,z)-G3(a,b-h,cc,z))/(2*h)
dC = (G3(a,b,cc+h,z)-G3(a,b,cc-h,z))/(2*h)
dZ = (G3(a,b,cc,z+h)-G3(a,b,cc,z-h))/(2*h)
# formula: D1=G(b,cc), D2=G(a,cc), D3=G(a,b)  (G2(a,b;z) etc.)
def G2(aa,bb,zz): return G_int([aa,bb],zz)
def G1(aa,zz): return mp.log(1-zz/aa)
D1, D2, D3 = G2(b,cc,z), G2(a,cc,z), G2(a,b,z)
# dG = D1 dlog((z-a)/(a-b)) + D2 dlog((a-b)/(b-cc)) + D3 dlog((b-cc)/cc)
# da component: D1*(-1/(z-a) - 1/(a-b)) + D2*(1/(a-b))
pred_dA = D1*(-1/(z-a) - 1/(a-b)) + D2*(1/(a-b))
pred_dB = D1*(1/(a-b)) + D2*(-1/(a-b) - 1/(b-cc)) + D3*(1/(b-cc))
pred_dC = D2*(1/(b-cc)) + D3*(-1/(b-cc) - 1/cc)
pred_dZ = D1*(1/(z-a))
print("dA:", mp.nstr(dA,20), " pred:", mp.nstr(pred_dA,20))
print("dB:", mp.nstr(dB,20), " pred:", mp.nstr(pred_dB,20))
print("dC:", mp.nstr(dC,20), " pred:", mp.nstr(pred_dC,20))
print("dZ:", mp.nstr(dZ,20), " pred:", mp.nstr(pred_dZ,20))
