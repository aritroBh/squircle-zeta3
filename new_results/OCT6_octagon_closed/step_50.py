import pickle, sympy as sp
r2 = sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
GENS = [l2,lu,z2,z3,L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h,L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h]
WT = {l2:1,lu:1,z2:2,z3:3,L2x1:2,L2x2:2,L2k:2,L2d:2,L2r2h:2,L2mr2h:2,L3x1:3,L3x2:3,L3k:3,L3d:3,L3r2h:3,L3mr2h:3}
def normQ2(q):
    num,den = sp.fraction(sp.together(q)); denc = den.xreplace({r2:-r2})
    return sp.collect(sp.expand(sp.expand(num*denc)/sp.expand(den*denc)), r2)
def monstr(mon):
    p=[]
    for g,e in zip(GENS,mon):
        if e: p.append(f"{g}"+(f"^{e}" if e>1 else ""))
    return '*'.join(p) if p else '1'
S = pickle.load(open('S_final.pkl','rb'))
def tab(expr):
    P = sp.Poly(sp.expand(expr), *GENS)
    rows=[]
    for mon,q in P.terms():
        w=sum(WT[g]*e for g,e in zip(GENS,mon))
        rows.append((w,monstr(mon),normQ2(q)))
    return sorted(rows,key=lambda r:(r[0],r[1]))
# P'''(1)/u - C4/u = 48*S_tot - (54-39r2)
e1 = sp.expand(48*S['Stot'] - (54-39*r2))
rows_diff = tab(e1)
print("P'''(1)/u - C4/u = 48*S_tot - (54-39r2):")
for w,m,q in rows_diff: print(f"  w{w} {m:14s} {q}")
pickle.dump(rows_diff, open('table_diffC4.pkl','wb'))
# numeric check
import mpmath as mp
mp.mp.dps=100
vals = pickle.load(open('basis745.pkl','rb'))
vmap = {l2:vals['l2'],lu:vals['lu'],z2:vals['z2'],z3:vals['z3'],
 L2x1:vals['Li2_x1'],L2x2:vals['Li2_x2'],L2k:vals['Li2_k1'],L2d:vals['Li2_d'],L2r2h:vals['Li2_r2h'],L2mr2h:vals['Li2_mr2h'],
 L3x1:vals['Li3_x1'],L3x2:vals['Li3_x2'],L3k:vals['Li3_k1'],L3d:vals['Li3_d'],L3r2h:vals['Li3_r2h'],L3mr2h:vals['Li3_mr2h']}
v = mp.mpf(sp.N(e1.subs(vmap, simultaneous=True),95))
P3q = mp.mpf(pickle.load(open('P3_exact.pkl','rb'))['P3_exact']); um=mp.sqrt(2+mp.sqrt(2))
print("numeric check vs P3/u-C4/u:", mp.nstr(abs(v-(P3q/um-(54-39*mp.sqrt(2)))),3))
