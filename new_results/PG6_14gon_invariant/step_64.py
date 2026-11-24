import pickle, sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
g = S.init()
S14 = pickle.load(open('S14.pkl','rb'))
# extract lone-atom terms: monomials with exactly one atom symbol, no letters/z
atom_terms = []  # (atom_name, coef_FT)
for mon,c in S14.items():
    d=dict(mon)
    syms = [s for s in d if s not in ('l2','la','lu','z2','z3')]
    others = [s for s in d if s in ('l2','la','lu','z2','z3')]
    if len(syms)==1 and not others:
        atom_terms.append((syms[0], c))
print("lone-atom terms:", len(atom_terms))
from collections import Counter
print(Counter(s for s,_ in atom_terms))
# map atom name -> (symbol, scalar)
def atom_sym(name):
    if name.startswith('L3x'): i=int(name[3:]); return S.sym_L3x(i)
    if name.startswith('Ax'): i=int(name[2:]); return S.sym_Ax(i)
    if name.startswith('Lam'):
        i,j = name[4],name[6]; return S.sym_Lam(int(i),int(j))
    if name.startswith('T3p'):
        i,j = name[4],name[6]; return S.sym_T3p(int(i),int(j))
    if name.startswith('T3_'):
        # T3_i_jk
        i=int(name[3]); j=int(name[5]); k=int(name[6]); return S.sym_T3a(i,j,k)
    raise ValueError(name)
# build sigma(S)_indecomp = sum coef * scalar * symbol
sigma={}
for name,c in atom_terms:
    sym,scalar = atom_sym(name)
    tot = E.fmul(c, scalar)   # FT coefficient
    for key,cf in sym.items():
        # store as (key) -> FT
        sigma[key] = E.fadd(sigma.get(key,E.ZERO), E.fs(tot, cf))
print("sigma(S)_indecomp terms:", len(sigma))
pickle.dump(sigma, open('sigma14.pkl','wb'))
# sanity: print a few atom symbol sizes
for nm in ['L3x1','Ax1','Lam1_2','T3p_1_2','T3_1_25']:
    sym,sc = atom_sym(nm)
    print(nm, "symbol terms:", len(sym), " scalar eval:", E.fev(sc,20))
open('progress.log','a').write('PG6 s_B8: computed 44 atom symbols + sigma(S)_indecomp (%d terms)\n'%len(sigma))
