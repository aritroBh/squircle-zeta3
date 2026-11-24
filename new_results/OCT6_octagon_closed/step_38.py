import pickle, sympy as sp
import mixed_engine as ME
S = pickle.load(open('S_assembly.pkl','rb'))
RED = pickle.load(open('reductions.pkl','rb'))
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
# map ME symbols to local canonical symbols
symmap = {ME.l2:l2, ME.lu:lu, ME.z2:z2, ME.z3:z3,
 ME.L2_1:L2x1, ME.L3_1:L3x1, ME.L2_2:L2x2, ME.L3_2:L3x2,
 ME.L2k:L2k, ME.L3k:L3k, ME.L2d:L2d, ME.L3d:L3d, ME.Lam:sp.symbols('LamS'), ME.swap:sp.symbols('swapS')}
# reduction map: ME.AA_1 -> R1, ME.AA_2 -> R2, Lam -> R3, swap -> R4 (R's already in canonical syms)
redmap = {ME.AA_1: RED['AA_1'], ME.AA_2: RED['AA_2'], sp.symbols('LamS'): RED['Lam'], sp.symbols('swapS'): RED['swap']}
Stot = S['S_tot'].subs(symmap, simultaneous=True)
Stot = sp.expand(Stot.subs(redmap, simultaneous=True))
Smix = S['S_mix'].subs(symmap, simultaneous=True)
Smix = sp.expand(Smix.subs(redmap, simultaneous=True))
print("S_tot after reduction: terms =", len(sp.Add.make_args(Stot)))
print("S_mix after reduction: terms =", len(sp.Add.make_args(Smix)))
pickle.dump({'Stot':Stot,'Smix':Smix}, open('S_final.pkl','wb'))
# quick numeric confirmation
import mpmath as mp, pickle as pk
mp.mp.dps = 745
vals = pk.load(open('basis745.pkl','rb'))
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
um = mp.sqrt(2+mp.sqrt(2))
vt = 48*um*mp.mpf(sp.N(Stot.subs(vmap, simultaneous=True),740))
vm = um*mp.mpf(sp.N(Smix.subs(vmap, simultaneous=True),740))
P3q = mp.mpf(pk.load(open('P3_exact.pkl','rb'))['P3_exact'])
print("48u*Stot_reduced vs quad P3:", mp.nstr(abs(vt-P3q),3))
print("u*Smix_reduced =", mp.nstr(vm,30), "vs target 0.00169887193388033285629748...")
