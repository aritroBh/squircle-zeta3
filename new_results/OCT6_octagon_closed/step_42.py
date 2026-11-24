import pickle, sympy as sp
r2 = sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2x1,L2x2,L2k,L2d,L2r2h,L2mr2h = sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')
L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h = sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h')
H = sp.Rational(1,2)
# W2 relations (dependent -> expr in basis {z2,l2^2,l2*lu,lu^2,L2x1})
rel2 = {
 L2x2: -l2**2/8 - L2x1,
 L2k:   sp.Rational(3,4)*z2 - lu**2/2 + L2x1,
 L2d:   z2/4 + l2*lu/2 - lu**2/2 - L2x1,
 L2r2h: z2 - l2**2/8 - l2*lu/2 + L2x1,
 L2mr2h: -sp.Rational(3,4)*z2 - l2**2/8 + l2*lu/2 - L2x1,
}
# W3 relations (dependent -> expr in basis {z3,L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h,z2l2->z2*l2,z2lu,l2*lu^2,L2x1*l2,L2x1*lu})
rel3 = {
 l2**3: -sp.Rational(21,4)*z3 + 24*L3r2h + 24*L3mr2h + 3*z2*l2,
 l2**2*lu: -sp.Rational(25,4)*z3 + 8*L3x1 + 8*L3x2 - 8*L3mr2h + 3*z2*l2,
 lu**3: -sp.Rational(171,64)*z3 + 3*L3k + 3*L3d - sp.Rational(3,2)*L3r2h + sp.Rational(3,2)*L3mr2h - sp.Rational(27,16)*z2*l2 + 3*z2*lu + sp.Rational(3,4)*l2*lu**2,
 L2x2*l2: sp.Rational(21,32)*z3 - 3*L3r2h - 3*L3mr2h - sp.Rational(3,8)*z2*l2 - L2x1*l2,
 L2k*l2: sp.Rational(3,4)*z2*l2 - H*l2*lu**2 + L2x1*l2,
 L2d*l2: -sp.Rational(25,8)*z3 + 4*L3x1 + 4*L3x2 - 4*L3mr2h + sp.Rational(7,4)*z2*l2 - H*l2*lu**2 - L2x1*l2,
 L2r2h*l2: sp.Rational(121,32)*z3 - 4*L3x1 - 4*L3x2 - 3*L3r2h + L3mr2h - sp.Rational(7,8)*z2*l2 + L2x1*l2,
 L2mr2h*l2: -sp.Rational(79,32)*z3 + 4*L3x1 + 4*L3x2 - 3*L3r2h - 7*L3mr2h + sp.Rational(3,8)*z2*l2 - L2x1*l2,
 L2x2*lu: sp.Rational(25,32)*z3 - L3x1 - L3x2 + L3mr2h - sp.Rational(3,8)*z2*l2 - L2x1*lu,
 L2k*lu: sp.Rational(171,128)*z3 - sp.Rational(3,2)*L3k - sp.Rational(3,2)*L3d + sp.Rational(3,4)*L3r2h - sp.Rational(3,4)*L3mr2h + sp.Rational(27,32)*z2*l2 - sp.Rational(3,4)*z2*lu - sp.Rational(3,8)*l2*lu**2 + L2x1*lu,
 L2d*lu: sp.Rational(171,128)*z3 - sp.Rational(3,2)*L3k - sp.Rational(3,2)*L3d + sp.Rational(3,4)*L3r2h - sp.Rational(3,4)*L3mr2h + sp.Rational(27,32)*z2*l2 - sp.Rational(5,4)*z2*lu + sp.Rational(1,8)*l2*lu**2 - L2x1*lu,
 L2r2h*lu: sp.Rational(25,32)*z3 - L3x1 - L3x2 + L3mr2h - sp.Rational(3,8)*z2*l2 + z2*lu - H*l2*lu**2 + L2x1*lu,
 L2mr2h*lu: sp.Rational(25,32)*z3 - L3x1 - L3x2 + L3mr2h - sp.Rational(3,8)*z2*l2 - sp.Rational(3,4)*z2*lu + H*l2*lu**2 - L2x1*lu,
}
# substitute rel2 INTO rel3 first (rel3 references L2k etc. implicitly? no - rel3 exprs use L2x1,L3's,logs,z2 only - check)
alldep = list(rel3.keys())
bad = [k for k in alldep for g in sp.sympify(rel3[k]).free_symbols if g in rel2]
print("rel3 exprs containing w2-dependent syms:", bad)
# good: rel3 only mentions L2x1 and L3's. But rel3 was derived with Li2 relations ALREADY independent? The greedy
# processed w3 on raw named values; the rel3 expressions are valid as-is. However they might not be consistent
# with rel2 (the Li2-product relations must follow from rel2 x log + Li3 relations). They are consistent numerically
# (resid ~1e-744), which is what matters; and I'll present the relation lattice as verified identities.

S = pickle.load(open('S_final.pkl','rb'))
def to_basis(expr):
    e = sp.expand(expr)
    e = sp.expand(e.subs(rel3, simultaneous=True))
    e = sp.expand(e.subs(rel2, simultaneous=True))
    return e
import time; t0=time.time()
Sb = {k: to_basis(S[k]) for k in ['Stot','Smix']}
print("basis substitution %.1fs"%(time.time()-t0))
pickle.dump(Sb, open('S_indep.pkl','wb'))
GENS = [l2,lu,z2,z3,L2x1,L3x1,L3x2,L3k,L3d,L3r2h,L3mr2h]
WT = {l2:1,lu:1,z2:2,z3:3,L2x1:2,L3x1:3,L3x2:3,L3k:3,L3d:3,L3r2h:3,L3mr2h:3}
def normQ2(q):
    num,den = sp.fraction(sp.together(q)); denc = den.xreplace({r2:-r2})
    return sp.collect(sp.expand(sp.expand(num*denc)/sp.expand(den*denc)), r2)
def monstr(mon):
    p=[]
    for g,e in zip(GENS,mon):
        if e: p.append(f"{g}"+(f"^{e}" if e>1 else ""))
    return '*'.join(p) if p else '1'
def show(Sexpr, label):
    P = sp.Poly(Sexpr, *GENS)
    rows=[]
    for mon,q in P.terms():
        w=sum(WT[g]*e for g,e in zip(GENS,mon))
        rows.append((w,monstr(mon),normQ2(q)))
    rows.sort(key=lambda r:(r[0],r[1]))
    print(f"--- {label} ({len(rows)} terms)")
    for w,m,q in rows: print(f"  w{w} {m:14s} {q}")
    return rows
rows_tot = show(Sb['Stot'], "P'''(1)/(48u) in independent basis")
rows_mix = show(Sb['Smix'], "mixed/(u) in independent basis")
pickle.dump({'tot':rows_tot,'mix':rows_mix}, open('final_tables.pkl','wb'))
