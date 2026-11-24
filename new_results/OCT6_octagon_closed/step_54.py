import pickle, sympy as sp, mpmath as mp
r2 = sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
sym_of = {'l2':l2,'lu':lu,'z2':z2,'z3':z3,'L2x1':L2x1,'L2x2':L2x2,'L2k':L2k,'L2d':L2d,'L2r2h':L2r2h,'L2mr2h':L2mr2h,
 'L3x1':L3x1,'L3x2':L3x2,'L3k':L3k,'L3d':L3d,'L3r2h':L3r2h,'L3mr2h':L3mr2h,'1':sp.Integer(1)}
def parse_monstr(m):
    if m=='1': return sp.Integer(1)
    e = sp.Integer(1)
    for f in m.split('*'):
        if '^' in f:
            b,pw = f.split('^'); e *= sym_of[b]**int(pw)
        else: e *= sym_of[f]
    return e
rows = pickle.load(open('table_Stot.pkl','rb'))
Stab = sum(q*parse_monstr(m) for w,m,q in rows)
mp.mp.dps = 745
vals = pickle.load(open('basis745.pkl','rb'))
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
um = mp.sqrt(2+mp.sqrt(2))
vt = 48*um*mp.mpf(sp.N(Stab.subs(vmap, simultaneous=True),740))
P3q = mp.mpf(pickle.load(open('P3_exact.pkl','rb'))['P3_exact'])
print("35-term table -> 48u*S vs quad P3:", mp.nstr(abs(vt-P3q),3))
# regenerate the md table lines programmatically for comparison
def fmt(w,m,q):
    qe = sp.collect(sp.expand(q), r2)
    return f"  w{w} {m}: {qe}"
print("\n".join(fmt(w,m,q) for w,m,q in rows))
