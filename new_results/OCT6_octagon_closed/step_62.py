import os, pickle, mpmath as mp, sympy as sp
# sanity: what's in oct5 (did my early steps overwrite OCT5's step_0X?) — originals are safe in OCT5's own dir; just report
o5 = sorted(f for f in os.listdir('oct5') if f.startswith('step_'))
print("oct5/ step files:", len(o5), o5[:4], "...")
# confirm local z3 arithmetic
print("48*(9/8 - 13*sqrt(2)/16) =", sp.expand(48*(sp.Rational(9,8) - 13*sp.sqrt(2)/16)))
print("(39*sqrt(2)-42)/4 =", sp.nsimplify((39*sp.sqrt(2)-42)/4), " ; 87/2-117*sqrt(2)/4 - (54-39*sqrt(2)) =", sp.expand(sp.Rational(87,2)-117*sp.sqrt(2)/4-(54-39*sp.sqrt(2))))
# final end-to-end from saved 35-term table
r2=sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
sym_of = {'l2':l2,'lu':lu,'z2':z2,'z3':z3,'L2x1':L2x1,'L2x2':L2x2,'L2k':L2k,'L2d':L2d,'L2r2h':L2r2h,'L2mr2h':L2mr2h,
 'L3x1':L3x1,'L3x2':L3x2,'L3k':L3k,'L3d':L3d,'L3r2h':L3r2h,'L3mr2h':L3mr2h,'1':sp.Integer(1)}
def parse_monstr(m):
    if m=='1': return sp.Integer(1)
    e=sp.Integer(1)
    for f in m.split('*'):
        if '^' in f: b,pw=f.split('^'); e*=sym_of[b]**int(pw)
        else: e*=sym_of[f]
    return e
mp.mp.dps=745
vals = pickle.load(open('basis745.pkl','rb'))
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
um=mp.sqrt(2+mp.sqrt(2))
rows = pickle.load(open('table_Stot.pkl','rb'))
v = 48*um*sum(mp.mpf(sp.N(q.subs(vmap),740))*mp.mpf(sp.N(parse_monstr(m).subs(vmap),740)) for w,m,q in rows)
P3q = mp.mpf(pickle.load(open('P3_exact.pkl','rb'))['P3_exact'])
print("\nFINAL end-to-end: 48u*(35-term named table) vs independent quad P3:", mp.nstr(abs(v-P3q),3))
print("P'''(1) =", mp.nstr(v,50))
with open('progress.log','a') as f:
    f.write("step_20: result.md finalized with analytic proofs of R1,R2. FINAL end-to-end table check: |diff|=%s. DONE.\n"%mp.nstr(abs(v-P3q),3))
print("logged")
