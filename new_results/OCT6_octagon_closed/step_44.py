import pickle, sympy as sp
import mixed_engine as ME
r2 = sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
S = pickle.load(open('S_assembly.pkl','rb'))
RED = pickle.load(open('reductions.pkl','rb'))
symmap = {ME.l2:l2, ME.lu:lu, ME.z2:z2, ME.z3:z3,
 ME.L2_1:L2x1, ME.L3_1:L3x1, ME.L2_2:L2x2, ME.L3_2:L3x2,
 ME.L2k:L2k, ME.L3k:L3k, ME.L2d:L2d, ME.L3d:L3d, ME.Lam:sp.symbols('LamS'), ME.swap:sp.symbols('swapS')}
redmap = {ME.AA_1: RED['AA_1'], ME.AA_2: RED['AA_2'], sp.symbols('LamS'): RED['Lam'], sp.symbols('swapS'): RED['swap']}
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
def tab(expr):
    P = sp.Poly(sp.expand(expr), *GENS)
    rows=[]
    for mon,q in P.terms():
        w=sum(WT[g]*e for g,e in zip(GENS,mon))
        rows.append((w,monstr(mon),normQ2(q)))
    return sorted(rows,key=lambda r:(r[0],r[1]))
out={}
for key in ['S_local','S_s1','S_s2']:
    e = S[key].subs(symmap, simultaneous=True)
    e = sp.expand(e.subs(redmap, simultaneous=True))
    rows = tab(e)
    out[key]=rows
    print(f"--- {key}  (x48u for P3)  z3-coeff:", [q for w,m,q in rows if m=='z3'])
    for w,m,q in rows: print(f"  w{w} {m:14s} {q}")
pickle.dump(out, open('subtables.pkl','wb'))
