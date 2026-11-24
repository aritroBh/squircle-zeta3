import pickle, sympy as sp, mpmath as mp
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
Sb = pickle.load(open('S_indep.pkl','rb'))
um = mp.sqrt(2+mp.sqrt(2))
vt = 48*um*mp.mpf(sp.N(Sb['Stot'].subs(vmap, simultaneous=True),740))
vm = um*mp.mpf(sp.N(Sb['Smix'].subs(vmap, simultaneous=True),740))
P3q = mp.mpf(pickle.load(open('P3_exact.pkl','rb'))['P3_exact'])
print("indep-basis 48u*Stot vs quad P3:", mp.nstr(abs(vt-P3q),3))
print("indep-basis u*Smix =", mp.nstr(vm,40))
# C4 difference
r2m = mp.sqrt(2)
C4 = (54-39*r2m)*um
print("\nC4 = (54-39r2)u =", mp.nstr(C4,30))
print("P3 - C4 =", mp.nstr(vt-C4,30), " (numeric)")
print("z3-coeff of P3/u (named 35-term form): 48*(29/32 - 39r2/64) =", mp.nstr(48*(mp.mpf(29)/32 - 39*r2m/64),20))
print("z3-coeff of P3/u (indep basis): 48*(-55891/512 + 39501r2/512) =", mp.nstr(48*(-mp.mpf(55891)/512 + 39501*r2m/512),20))
