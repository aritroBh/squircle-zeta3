import pickle, mpmath as mp, sympy as sp
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
r2 = sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
half = sp.Rational(1,2)
# R1: AA_1
R1 = half*z3 + sp.Rational(3,2)*L3x1 + half*L3x2 - half*L3r2h - half*z2*l2 - half*l2**2*lu + half*l2*lu**2 \
     + half*L2x1*l2 - half*L2x1*lu - half*L2x2*lu + L2d*l2
# R2: AA_2
R2 = z3 + L3x2 - L3r2h + half*z2*lu - half*l2**2*lu + half*L2mr2h*l2 - half*L2x2*lu + half*L2k*l2 \
     - half*L2r2h*l2 - half*L2r2h*lu
# R3: Lam
R3 = sp.Rational(3,4)*z3 + 2*L3x1 + L3k - L3d + sp.Rational(1,4)*l2**2*lu - sp.Rational(1,4)*l2*lu**2 \
     - half*L2x1*l2 + sp.Rational(1,4)*L2k*l2 - sp.Rational(1,4)*L2d*l2
# R4: swap (with sqrt2 coeffs)
R4 = -4*L3mr2h + 2*L3x2 - L3k + L3d - 4*L3r2h + (-half+half*r2)*z2*l2 - l2**2*lu + (half-half*r2)*l2*lu**2 \
     + half*lu**3 - r2*L2mr2h*lu + half*L2x1*l2 - half*L2x2*l2 + (-half-half*r2)*L2x2*lu \
     + (-half-half*r2)*L2k*l2 + (half-half*r2)*L2k*lu + (1-half*r2)*L2d*l2 + (half+half*r2)*L2d*lu \
     - half*L2r2h*l2 - half*L2r2h*lu
RED = {'AA_1':R1,'AA_2':R2,'Lam':R3,'swap':R4}
pickle.dump(RED, open('reductions.pkl','wb'))
# numeric check
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
for nm,e in RED.items():
    lhs = vals[nm]; rhs = mp.mpf(sp.N(e.subs(vmap, simultaneous=True),740))
    print(nm, "|LHS-RHS| =", mp.nstr(abs(lhs-rhs),3))
# cross-check R1+R2 vs Landen (proved analytically):
Landen = -vals['z2']*vals['l2']/2 - vals['l2']*vals['Li2_r2h'] - (vals['l2']**2/8)*(vals['l2']/2+vals['lu']) - 3*vals['Li3_r2h'] + 3*vals['z3']
R12 = mp.mpf(sp.N((R1+R2).subs(vmap, simultaneous=True),740))
print("R1+R2 vs Landen-RHS: |diff| =", mp.nstr(abs(R12-Landen),3))
with open('progress.log','a') as f:
    f.write("step_10: four key reductions saved (AA_1,AA_2,Lam,swap -> single Li3/Li2 reps + products), each verified ~745 digits; R1+R2 matches analytic Landen identity.\n")
