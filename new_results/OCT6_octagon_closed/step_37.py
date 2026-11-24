import pickle, mpmath as mp, sympy as sp
mp.mp.dps = 780
r2=mp.sqrt(2); k1=r2-1; x1=-k1
def Li21N(a,b,N):
    s=mp.mpf(0); ap=a*a; inner=b; bp=b
    for n in range(2,N+1):
        if n>2: ap*=a; bp*=b; inner+=bp/(n-1)
        s += ap*inner/n**2
    return s
Lam_best = Li21N(x1,-r2/2,2200)   # 0.414^2200 ~ 1e-847
Lam_chk  = Li21N(x1,-r2/2,1800)   # 0.414^1800 ~ 1e-696
print("|Lam2200 - Lam1800| =", mp.nstr(abs(Lam_best-Lam_chk),3), "(expect ~1e-696)")
vals = pickle.load(open('basis745.pkl','rb'))
print("|Lam_best - stored| =", mp.nstr(abs(Lam_best-vals['Lam']),3))
vals['Lam']=Lam_best
pickle.dump(vals, open('basis745.pkl','wb'))
# re-verify R3
RED = pickle.load(open('reductions.pkl','rb'))
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
for nm in ['AA_1','AA_2','Lam','swap']:
    rhs = mp.mpf(sp.N(RED[nm].subs(vmap, simultaneous=True),740))
    print(nm, "|LHS-RHS| =", mp.nstr(abs(vals[nm]-rhs),3))
