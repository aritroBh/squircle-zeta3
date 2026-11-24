import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
mp.mp.dps = 620
tl = pickle.load(open('dod5_tallies.pkl','rb'))
Tot3, TotP, TotI = tl['Tot3'], tl['TotP'], tl['TotI']
ZERO=(0,0)
print("=== Total Li3 coefficients over the 8-arg basis (in S) ===")
cnt = 0
for j in range(8):
    if Tot3[j] != ZERO:
        cnt += 1
        print("   Li3(%-12s ~%+.8f): (%s) + (%s) r3"%(str(BASIS3_ZS[j]), pfloat(BASIS3_ZS[j]), Tot3[j][0], Tot3[j][1]))
    else:
        print("   Li3(%-12s ~%+.8f): ZERO"%(str(BASIS3_ZS[j]), pfloat(BASIS3_ZS[j])))
print("nonzero basis args:", cnt, "of 8")
print("\n=== Total invisible coefficients (F-basis, in S) ===")
nms = ['zeta(3)','L(3,chi12)','zeta(2)l2','zeta(2)l3','zeta(2)lu']
for m in range(5):
    print("   %-12s: (%s) + (%s) r3"%(nms[m], TotI[m][0], TotI[m][1]))
# the local zeta(3) comparison
from fractions import Fraction as Fr
loc = (Fr(1,16), Fr(-3,32))
print("\nlocal part (2-3r3)/32 = (1/16) + (-3/32) r3;  total z3 coeff = (%s)+(%s)r3"%(TotI[0][0], TotI[0][1]))
nonloc = (TotI[0][0]-loc[0], TotI[0][1]-loc[1])
print("non-local z3 part:", nonloc)
# P'''(1) level: multiply by 72L
r3m = mp.sqrt(3); L = mp.sqrt(26-15*r3m)
def fval(p):
    return mp.mpf(p[0].numerator)/p[0].denominator + r3m*mp.mpf(p[1].numerator)/p[1].denominator
print("\n72L*(total z3 coeff) =", mp.nstr(72*L*fval(TotI[0]), 25))
print("C6 = (9/4)(2-3r3)L   =", mp.nstr(Fr(9,4)*(2-3*r3m)*L, 25))
print("72L*(L3chi coeff)    =", mp.nstr(72*L*fval(TotI[1]), 25))
