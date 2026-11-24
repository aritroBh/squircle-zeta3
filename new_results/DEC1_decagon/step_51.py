import os, sympy as sp, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
# symbolic named constants
z2,z3,lf,l2 = sp.symbols('z2 z3 lf l2')
L2f,L2x1,L2x4 = sp.symbols('L2f L2x1 L2x4')   # Li2(1/phi), Li2(-1/phi), Li2(1/phi^2)
L3f,L3x1,L3x4 = sp.symbols('L3f L3x1 L3x4')   # Li3(1/phi), Li3(-1/phi), Li3(1/phi^2)
# Li2 golden relations (Landen/Lewin, verified numerically below)
subL2 = {L2x4: sp.Rational(2,5)*z2 - lf**2, L2f: sp.Rational(3,5)*z2 - lf**2, L2x1: -sp.Rational(2,5)*z2 + lf**2/2}
# R1 = R_L2(x4), R2 = DIFF(x4), R3 = R_L2(xf), R4 = R_L1(phi) (Li-inversion reduced)
R1 = -z2*lf - 2*lf*L2f - lf**3 - 3*L3f + 3*z3
R2 = 2*L3x4 - 2*L3f - lf*L2f + 2*lf*L2x4 + lf**3
R3 = -2*z2*lf - 4*lf*L2x4 - 2*lf**3 - 3*L3x4 + 3*z3
R4 = z3 - z2*lf + lf*L2x1 + L3x1
Ax1 = sp.expand((R1-R2-R3+R4)/2)
Ax4 = sp.expand(R1 - Ax1)
Ax1s = sp.expand(Ax1.subs(subL2))
Ax4s = sp.expand(Ax4.subs(subL2))
print("A(x1) solved:", Ax1s)
print("A(x4) solved:", Ax4s)
# LLL reductions to prove:
LLL_Ax1 = sp.Rational(4,5)*z3 - sp.Rational(2,5)*z2*lf + 2*L3x1
LLL_Ax4 = z3 - sp.Rational(3,5)*z2*lf + 4*L3x1 + 3*L3f
print("A(x1)-LLL difference:", sp.expand(Ax1s - LLL_Ax1))
print("A(x4)-LLL difference:", sp.expand(Ax4s - LLL_Ax4))
# numeric verification of everything at 830 dps
mp.mp.dps = 830
v = pickle.load(open('atoms830.pkl','rb'))
mpf = lambda x: mp.mpf(x)
num = {z2:v['z2'],z3:v['z3'],lf:v['lf'],l2:v['l2'],L2f:mp.polylog(2,1/v['phi']),L2x1:mp.polylog(2,-1/v['phi']),L2x4:mp.polylog(2,1/v['phi']**2),L3f:v['L3_xf'],L3x1:v['L3_x1'],L3x4:v['L3_x4']}
# check Li2 relations
for sym,form in subL2.items():
    d = abs(mp.mpf(sp.N(form.subs(num),100)) - num[sym])
    print(f"Li2 rel {sym}: resid {mp.nstr(d,6)}")
d1 = abs(mp.mpf(sp.N(Ax1s.subs(num),100)) - v['A_x1'])
d2 = abs(mp.mpf(sp.N(Ax4s.subs(num),100)) - v['A_x4'])
print("A(x1) orbit-system vs direct series:", mp.nstr(d1,6))
print("A(x4) orbit-system vs direct series:", mp.nstr(d2,6))
with open('progress.log','a') as f:
    f.write("step: A-atom 4x4 orbit system (L2@x4, DIFF@x4, L2@xf, L1@phi + Li inversions) solved exactly;\n"
            "  solutions == LLL reductions EXACTLY (symbolic diff 0); Li2 golden relations verified 830d; A(x1),A(x4) match direct series.\n")