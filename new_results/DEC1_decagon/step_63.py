import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
S3 = pickle.load(open('S_reduced.pkl','rb'))['S3']
lf,l2,z2,z3 = E.lf,E.l2,E.z2,E.z3
L3_1,L3f = E.L3_1,E.L3f
vars6 = [l2,lf,z2,z3,L3_1,L3f]
# collect coefficients of monomials in vars6
poly = sp.Poly(sp.expand(S3), *vars6)
table = {}
for monom, coeff in poly.terms():
    if coeff==0: continue
    table[monom] = sp.expand(coeff)
print("distinct named-constant monomials:", len(table))
# pretty print
names = {l2:'l2',lf:'lf',z2:'z2',z3:'z3',L3_1:'L3(-1/phi)',L3f:'L3(1/phi)'}
def mname(m):
    parts=[]
    for v,e in zip(vars6,m):
        if e==1: parts.append(names[v])
        elif e>1: parts.append(f"{names[v]}^{e}")
    return '*'.join(parts) if parts else '1'
tot = sp.Integer(0)
for m in sorted(table, key=lambda mm:(sum(mm),mm)):
    c = table[m]
    tot += c
    print(f"{mname(m):24s} : {c}")
print("rational part:", tot)
pickle.dump({'table':{mname(m):table[m] for m in table}}, open('final_table.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: final table: {len(table)} named monomials in basis [l2,lf,z2,z3,L3(-1/phi),L3(1/phi)]; saved final_table.pkl\n")