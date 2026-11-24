import pickle, sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
g = S.init()
S14 = pickle.load(open('S14.pkl','rb'))
def is_w3_atom(nm):
    return nm.startswith('L3x') or nm.startswith('Ax') or nm.startswith('Lam') or nm.startswith('T3p') or nm.startswith('T3_')
atom_terms=[]
for mon,c in S14.items():
    syms=[s for s in dict(mon) if s not in ('l2','la','lu','z2','z3')]
    others=[s for s in dict(mon) if s in ('l2','la','lu','z2','z3')]
    if len(syms)==1 and not others and is_w3_atom(syms[0]):
        atom_terms.append((syms[0],c))
def atom_sym(name):
    if name.startswith('L3x'): return S.sym_L3x(int(name[3:]))
    if name.startswith('Ax'): return S.sym_Ax(int(name[2:]))
    if name.startswith('Lam'): return S.sym_Lam(int(name[3]),int(name[5]))
    if name.startswith('T3p'): return S.sym_T3p(int(name[4]),int(name[6]))
    if name.startswith('T3_'): return S.sym_T3a(int(name[3]),int(name[5]),int(name[6]))
    raise ValueError(name)
# test each atom type once
for nm in ['L3x1','Ax1','Lam1_2','T3p_1_2','T3_1_25']:
    sym,sc=atom_sym(nm); print(nm,"sym terms",len(sym))
sigma={}
for name,c in atom_terms:
    sym,scalar = atom_sym(name)
    tot=E.fmul(c,scalar)
    for key,cf in sym.items():
        sigma[key]=E.fadd(sigma.get(key,E.ZERO), E.fs(tot,cf))
print("weight-3 lone atoms:", len(atom_terms), " sigma terms:", len(sigma))
pickle.dump(sigma, open('sigma14.pkl','wb'))
open('progress.log','a').write('PG6 s_B10: sigma(S)_w3 = %d terms from %d weight-3 atoms\n'%(len(sigma),len(atom_terms)))
