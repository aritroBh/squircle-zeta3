import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
r5 = sp.sqrt(5)
tab = pickle.load(open('final_table_pq.pkl','rb'))['table_pq']
# high-precision eval of the 15-term table at 700 digits
mp.mp.dps = 750
v = pickle.load(open('atoms830.pkl','rb'))
lf,l2,z2,z3 = v['lf'],v['l2'],v['z2'],v['z3']
vals = {'l2':l2,'lf':lf,'z2':z2,'z3':z3,'L3(-1/phi)':v['L3_x1'],'L3(1/phi)':v['L3_xf']}
S = mp.mpf(0)
for mname,(p,q) in tab.items():
    mval = mp.mpf(1)
    if mname!='1':
        for part in mname.split('*'):
            if '^' in part:
                nm,e = part.split('^'); mval *= vals[nm]**int(e)
            else:
                mval *= vals[part]
    S += (mp.mpf(p)+mp.mpf(q)*mp.sqrt(5))*mval
s36 = mp.sqrt(10-2*mp.sqrt(5))/4; Lm = 2*s36/(mp.sqrt(5)*v['phi']**2)
P3_final = 60*Lm*S
gold = pickle.load(open('P3_660.pkl','rb'))['P3']
print("|P3 (15-term exact) - 640d gold| =", mp.nstr(abs(P3_final-gold),8))
# clean summary quantities
z3t = sp.Rational(33,64)-sp.Rational(19,80)*r5
print("[z3] S3 =", z3t)
T_over_L = sp.simplify(60*z3t)
print("total z3 / L =", T_over_L)
C5_over_L = sp.Rational(75,8)-sp.Rational(15,2)*r5
print("C_5 / L =", C5_over_L)
Ne_over_L = sp.simplify((T_over_L - C5_over_L)/10)
print("N_e / L =", Ne_over_L, "= 3(115-36 sqrt5)/160 =", sp.Rational(3,160)*(115-36*r5))
print("L =", sp.nsimplify(Lm) if False else "2 sin36/(sqrt5 phi^2)")
print("L^2 = (25-11 sqrt5)/10")
with open('progress.log','a') as f:
    f.write(f"step: 15-term table at 750dps vs 640d gold: {mp.nstr(abs(P3_final-gold),4)}. z3: [S3]={z3t}, T/L={T_over_L}, C5/L={C5_over_L}, Ne/L={Ne_over_L}\n")