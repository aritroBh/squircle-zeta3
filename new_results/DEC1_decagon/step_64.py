import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
S3 = pickle.load(open('S_reduced.pkl','rb'))['S3']
r5 = sp.sqrt(5)
lf,l2,z2,z3 = E.lf,E.l2,E.z2,E.z3
L3_1,L3f = E.L3_1,E.L3f
vars6 = [l2,lf,z2,z3,L3_1,L3f]
poly = sp.Poly(sp.expand(S3), *vars6)
def pq(c):
    c2 = sp.expand(c)
    p = sp.simplify((c2 + c2.subs(r5,-r5))/2)
    q = sp.simplify((c2 - c2.subs(r5,-r5))/(2*r5))
    return p,q
table = {}
for monom, coeff in poly.terms():
    if coeff==0: continue
    table[monom] = pq(coeff)
names = {l2:'l2',lf:'lf',z2:'z2',z3:'z3',L3_1:'L3(-1/phi)',L3f:'L3(1/phi)'}
def mname(m):
    parts=[]
    for v,e in zip(vars6,m):
        if e==1: parts.append(names[v])
        elif e>1: parts.append(f"{names[v]}^{e}")
    return '*'.join(parts) if parts else '1'
print("distinct named-constant monomials:", len(table))
tot = sp.Integer(0)
for m in sorted(table, key=lambda mm:(sum(mm),mm)):
    p,q = table[m]
    tot += p
    print(f"{mname(m):16s} : ({p}) + ({q}) sqrt5")
print("sum of rational parts (should be 0):", sp.simplify(tot))
pickle.dump({'table_pq':{mname(m):table[m] for m in table}}, open('final_table_pq.pkl','wb'))
# numeric check of the table
mp.mp.dps=110
v = pickle.load(open('atoms830.pkl','rb'))
num = {E.l2:v['l2'], E.lf:v['lf'], E.z2:v['z2'], E.z3:v['z3'], E.L3_1:v['L3_x1'], E.L3f:v['L3_xf']}
S3n = mp.mpf(sp.N(S3.subs(num, simultaneous=True), 110))
S_tab = mp.mpf(0)
for m,(p,q) in table.items():
    val = mp.mpf(1)
    for vv,e in zip(vars6,m):
        if e: val *= num[vv]**e
    S_tab += (mp.mpf(p)+mp.mpf(q)*mp.sqrt(5))*val
print("table vs S3:", mp.nstr(abs(S_tab-S3n),8))
with open('progress.log','a') as f:
    f.write(f"step: canonical table: {len(table)} monomials, p+q sqrt5 coefficients extracted; matches S3.\n")