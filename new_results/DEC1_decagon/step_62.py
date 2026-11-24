import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
S = pickle.load(open('S_assembly.pkl','rb'))['S']
lf,l2,z2,z3 = E.lf,E.l2,E.z2,E.z3
L3_1,L3f,L3_4 = E.L3_1,E.L3f,E.L3_4
L2_1,L2_4,L2f = E.L2_1,E.L2_4,E.L2f
# reductions (PROVED for A-atoms; LLL+700d cross-validated for Lam/swap)
red = {
 E.AA_1: sp.Rational(4,5)*z3 - sp.Rational(2,5)*z2*lf + 2*L3_1,
 E.AA_4: z3 - sp.Rational(3,5)*z2*lf + 4*L3_1 + 3*L3f,
 E.Lam:  sp.Rational(2,5)*z3 + sp.Rational(4,5)*z2*lf - lf**3 + 3*L3_1 + L3f,
 E.swap: sp.Rational(3,5)*z3 - sp.Rational(7,5)*z2*lf + 2*lf**3 + L3_1 + L3f,
}
S1 = sp.expand(S.subs(red, simultaneous=True))
print("terms after atom reduction:", len(sp.Add.make_args(S1)))
# Li2 golden reductions (classical, verified 830d)
red2 = {L2_1: -sp.Rational(2,5)*z2 + lf**2/2, L2_4: sp.Rational(2,5)*z2 - lf**2, L2f: sp.Rational(3,5)*z2 - lf**2}
S2 = sp.expand(S1.subs(red2, simultaneous=True))
print("terms after Li2 reduction:", len(sp.Add.make_args(S2)))
# Li3(1/phi^2) via R_a (keeps L3x1, L3xf):  L3_4 = 4 L3_1 + 4 L3f
S3 = sp.expand(S2.subs({L3_4: 4*L3_1+4*L3f}, simultaneous=True))
print("terms after R_a:", len(sp.Add.make_args(S3)))
pickle.dump({'S1':S1,'S2':S2,'S3':S3}, open('S_reduced.pkl','wb'))
# numeric verify S3 at high precision
mp.mp.dps = 110
v = pickle.load(open('atoms830.pkl','rb'))
num = {E.l2:v['l2'], E.lf:v['lf'], E.z2:v['z2'], E.z3:v['z3'], E.L3_1:v['L3_x1'], E.L3f:v['L3_xf']}
S3n = mp.mpf(sp.N(S3.subs(num, simultaneous=True), 110))
s36 = mp.sqrt(10-2*mp.sqrt(5))/4; Lm = 2*s36/(mp.sqrt(5)*v['phi']**2)
gold = pickle.load(open('P3_660.pkl','rb'))['P3']
print("|60 L S3 - gold| =", mp.nstr(abs(60*Lm*S3n-gold),6))
with open('progress.log','a') as f:
    f.write(f"step: reductions substituted: S1 ({len(sp.Add.make_args(S1))} terms), S2 ({len(sp.Add.make_args(S2))}), S3 ({len(sp.Add.make_args(S3))}); |60 L S3 - gold|={mp.nstr(abs(60*Lm*S3n-gold),4)}\n")