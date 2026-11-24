import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
import dec_engine as E
r5 = sp.sqrt(5)
d = pickle.load(open('S_assembly_dec2.pkl','rb'))
S3 = d['S3']
vars6 = [E.l2, E.lf, E.z2, E.z3, E.L3_1, E.L3f]
poly = sp.Poly(sp.expand(S3), *vars6)
def pq(c):
    c2 = sp.together(c)
    cc = c2.subs(r5,-r5)
    p = sp.radsimp(sp.expand((c2+cc)/2))
    q = sp.radsimp(sp.expand((c2-cc)/(2*r5)))
    return p,q
names = {E.l2:'l2', E.lf:'lf', E.z2:'z2', E.z3:'z3', E.L3_1:'L3(-1/phi)', E.L3f:'L3(1/phi)'}
def mname(m):
    parts=[]
    for v,e in zip(vars6,m):
        if e==1: parts.append(names[v])
        elif e>1: parts.append(f"{names[v]}^{e}")
    return '*'.join(parts) if parts else '1'
table = {}
for monom, coeff in poly.terms():
    if sp.expand(coeff)==0: continue
    table[mname(monom)] = pq(coeff)
print("my table terms:", len(table))
ft = pickle.load(open('<project>/research/08_theory/DEC1_decagon/final_table_pq.pkl','rb'))['table_pq']
print("DEC1 table terms:", len(ft))
mism = 0
for k in ft:
    if k not in table:
        print("MISSING", k); mism += 1; continue
    p1 = sp.simplify(sp.sympify(table[k][0]) - sp.sympify(ft[k][0]))
    p2 = sp.simplify(sp.sympify(table[k][1]) - sp.sympify(ft[k][1]))
    if p1 != 0 or p2 != 0:
        mism += 1
        print("MISMATCH", k, "mine:", table[k], " dec1:", ft[k])
extra = set(table)-set(ft)
print("keys only in mine:", extra)
print("mismatches:", mism)
print("rational term '1' in my table:", table.get('1', "ABSENT (correct)"))
for k in sorted(table, key=lambda s:(s.count('*'),s)):
    p,q = table[k]
    print(f"  {k:16s} : ({p}) + ({q}) sqrt5")
pickle.dump(table, open('final_table_dec2.pkl','wb'))
