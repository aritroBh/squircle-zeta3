import os, pickle, time
import sympy as sp
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'E5_polys.pkl'), 'rb') as f: E5_polys = pickle.load(f)
with open(os.path.join(k1f, 'jtab5.pkl'), 'rb') as f: jtab5 = pickle.load(f)
with open(os.path.join(k1f, 'jloc_sym.pkl'), 'rb') as f: jloc = pickle.load(f)
E5 = E5_polys[5]
x = sp.Symbol('x')
atom_names = ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5',
              'H00011','H00101','H00111','H01011','H01111']
atoms = {nm: sp.Symbol(nm) for nm in atom_names}
z2,z3,z5 = sp.symbols('z2 z3 z5')
ns_nl = dict(atoms); ns_nl['x']=x; ns_nl['z2']=z2; ns_nl['z3']=z3; ns_nl['z5']=z5
J5 = {k: sp.sympify(s, locals=ns_nl) for k,s in jtab5.items()}
DEN = sp.Symbol('DEN'); HH = {i: sp.Symbol(f'HH{i}') for i in range(1,6)}
ns_loc = {'DEN':DEN,'z2':z2,'z3':z3,'z5':z5}
for i in range(1,6): ns_loc[f'H{i}']=HH[i]
Jloc = {k: sp.sympify(s, locals=ns_loc) for k,s in jloc.items()}
def localJ(n,a,b):
    subs = {DEN: sp.Rational(1,n+1)}
    for k in range(1,6): subs[HH[k]] = sp.harmonic(n+1,k)
    return sp.expand(Jloc[(a,b)].subs(subs))
JLOC_n = {}
for (n,a,b,c) in [k for k in E5 if k[3]==0]:
    JLOC_n[(n,a,b)] = localJ(n,a,b)

# Save these building blocks for reuse
import pickle as pk
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/build_blocks.pkl','wb') as f:
    pk.dump({'J5':J5,'JLOC_n':JLOC_n,'atoms':atoms,'z2':z2,'z3':z3,'z5':z5,'x':x}, f)
print("Saved build_blocks.pkl")
print("num J5:", len(J5), " num JLOC_n:", len(JLOC_n))
