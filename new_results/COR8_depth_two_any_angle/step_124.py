
import sys, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols_2kept.pkl','rb'))
OM, DT, R = D['OM'], D['DT'], D['R']
dR = delta22(R, 4)
print("delta22(Omega-DT) with [2] kept:")
names = ['2','c','m','p']
for k, poly in sorted(dR.items()):
    (i1,j1),(i2,j2) = k
    print(f"  ({names[i1]}^{names[j1]})^({names[i2]}^{names[j2]}):", {e: str(v) for e, v in sorted(poly.items())})
print("delta22(DT):")
dDT = delta22(DT, 4)
for k, poly in sorted(dDT.items()):
    (i1,j1),(i2,j2) = k
    print(f"  ({names[i1]}^{names[j1]})^({names[i2]}^{names[j2]}):", {e: str(v) for e, v in sorted(poly.items())})
