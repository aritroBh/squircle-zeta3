import pickle, sympy as sp
S = pickle.load(open('S_final.pkl','rb'))
r2 = sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
L2S = list(sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h'))
L3S = list(sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h'))
GENS = [l2,lu,z2,z3]+L2S+L3S
WT = {l2:1,lu:1,z2:2,z3:3, **{g:2 for g in L2S}, **{g:3 for g in L3S}}
def normQ2(q):
    num,den = sp.fraction(sp.together(q))
    denc = den.xreplace({r2:-r2})
    e = sp.expand(num*denc)/sp.expand(den*denc)
    return sp.collect(sp.expand(e), r2)
def monstr(mon):
    parts=[]
    for g,e in zip(GENS,mon):
        if e: parts.append(f"{g}"+(f"^{e}" if e>1 else ""))
    return '*'.join(parts) if parts else '1'
def table(Sexpr):
    P = sp.Poly(sp.expand(Sexpr), *GENS)
    rows=[]
    for mon,q in P.terms():
        w = sum(WT[g]*e for g,e in zip(GENS,mon))
        rows.append((w, monstr(mon), normQ2(q)))
    return sorted(rows, key=lambda r:(r[0],r[1]))
for name in ['Stot','Smix']:
    rows = table(S[name])
    print("="*100)
    print(f"{name}:  P'''/{'48u' if name=='Stot' else 'u'} = Σ q·m   ({len(rows)} terms)")
    for w,m,q in rows:
        print(f"  w{w}  {m:22s} {q}")
    pickle.dump(rows, open(f'table_{name}.pkl','wb'))
