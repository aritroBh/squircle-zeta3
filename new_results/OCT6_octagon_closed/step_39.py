import pickle, sympy as sp
Stot = pickle.load(open('S_final.pkl','rb'))['Stot']
r2 = sp.sqrt(2)
l2,lu,z2,z3 = sp.symbols('l2 lu z2 z3')
GENS = [l2,lu,z2,z3] + list(sp.symbols('L2x1 L2x2 L2k L2d L2r2h L2mr2h')) + list(sp.symbols('L3x1 L3x2 L3k L3d L3r2h L3mr2h'))
P = sp.Poly(Stot, *GENS)
rows = []
for mon, q in P.terms():
    w = sum(mon)  # total degree = weight
    rows.append((w, mon, q))
rows.sort()
def monstr(mon):
    parts=[]
    for g,e in zip(GENS,mon):
        if e: parts.append(f"{g}" + (f"^{e}" if e>1 else ""))
    return '*'.join(parts) if parts else '1'
print("P'''(1)/(48 u) = sum of q_m * m   [weight, monomial m, q in Q(sqrt2)]")
for w,mon,q in rows:
    print(f"  w{w}  {monstr(mon):28s}  {sp.nsimplify(q)}")
print("\nnumber of nonzero monomials:", len(rows))
