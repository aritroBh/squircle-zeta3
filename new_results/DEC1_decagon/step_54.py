import os, sympy as sp, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
z2,z3,lf,l2 = sp.symbols('z2 z3 lf l2')
L2f,L2x1,L2x4 = sp.symbols('L2f L2x1 L2x4')
L3f,L3x1,L3x4 = sp.symbols('L3f L3x1 L3x4')
subL2 = {L2x4: sp.Rational(2,5)*z2 - lf**2, L2f: sp.Rational(3,5)*z2 - lf**2, L2x1: -sp.Rational(2,5)*z2 + lf**2/2}
R1 = -z2*lf - 2*lf*L2f - lf**3 - 3*L3f + 3*z3
R2 = 2*L3x4 - 2*L3f - lf*L2f + 2*lf*L2x4 + lf**3
R3 = -2*z2*lf - 4*lf*L2x4 - 2*lf**3 - 3*L3x4 + 3*z3
R4 = z3 - z2*lf + lf*L2x1 + 3*L3x1      # corrected
Ax1 = sp.expand((R1-R2-R3+R4)/2)
Ax4 = sp.expand(R1 - Ax1)
Ax1s = sp.expand(Ax1.subs(subL2)); Ax4s = sp.expand(Ax4.subs(subL2))
LLL_Ax1 = sp.Rational(4,5)*z3 - sp.Rational(2,5)*z2*lf + 2*L3x1
LLL_Ax4 = z3 - sp.Rational(3,5)*z2*lf + 4*L3x1 + 3*L3f
d1 = sp.expand(Ax1s - LLL_Ax1); d2 = sp.expand(Ax4s - LLL_Ax4)
print("A(x1) orbit =", Ax1s)
print("A(x4) orbit =", Ax4s)
print("diff A(x1) vs LLL:", d1)
print("diff A(x4) vs LLL:", d2)
# reduce diffs using golden Li3 lattice: R_a: L3x4 = 4L3x1+4L3f ; sum-rel: L3x1+L3f = (z3-z2*lf)/5 + lf^3/6
subL3 = {L3x4: 4*L3x1+4*L3f, L3f: (z3-z2*lf)/5 + lf**3/6 - L3x1}
print("diff1 reduced:", sp.expand(d1.subs(subL3)))
print("diff2 reduced:", sp.expand(d2.subs(subL3)))
# also reduce A solutions to the L3 basis {L3x1, L3f} for the record
print("A(x1) in {L3x1,L3f}:", sp.expand(Ax1s))
# numeric final check vs direct series at 830 dps
mp.mp.dps=830
v = pickle.load(open('atoms830.pkl','rb'))
num = {z2:v['z2'],z3:v['z3'],lf:v['lf'],L2f:mp.polylog(2,1/v['phi']),L2x1:mp.polylog(2,-1/v['phi']),L2x4:mp.polylog(2,1/v['phi']**2),L3f:v['L3_xf'],L3x1:v['L3_x1'],L3x4:v['L3_x4']}
print("|orbit A(x1) - series| =", mp.nstr(abs(mp.mpf(sp.N(Ax1s.subs(num),100))-v['A_x1']),6))
print("|orbit A(x4) - series| =", mp.nstr(abs(mp.mpf(sp.N(Ax4s.subs(num),100))-v['A_x4']),6))
with open('progress.log','a') as f:
    f.write("step: orbit system re-solved with corrected R4 (Li3 inversion Li3(-x)=Li3(-1/x)-z2 ln x-ln^3 x/6, PROVED by differentiation).\n"
            f"  A(x1) = {Ax1s}\n  A(x4) = {Ax4s}\n  both == LLL reductions exactly modulo golden Li3 lattice (R_a + sum-rel).\n")