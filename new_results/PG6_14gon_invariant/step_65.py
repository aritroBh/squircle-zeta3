import pickle, sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
g = S.init()
S14 = pickle.load(open('S14.pkl','rb'))
W3 = ('L3x','Ax','Lam','T3p','T3_')
def is_w3_atom(nm):
    return nm.startswith('L3x') or nm.startswith('Ax') or nm.startswith('Lam') or nm.startswith('T3p') or nm.startswith('T3_')
# weight-3 lone atoms
atom_terms=[]
for mon,c in S14.items():
    syms=[s for s in dict(mon) if s not in ('l2','la','lu','z2','z3')]
    others=[s for s in dict(mon) if s in ('l2','la','lu','z2','z3')]
    if len(syms)==1 and not others and is_w3_atom(syms[0]):
        atom_terms.append((syms[0],c))
print("weight-3 lone atoms:", len(atom_terms))
def atom_sym(name):
    if name.startswith('L3x'): return S.sym_L3x(int(name[3:]))
    if name.startswith('Ax'): return S.sym_Ax(int(name[2:]))
    if name.startswith('Lam'): return S.sym_Lam(int(name[4]),int(name[6]))
    if name.startswith('T3p'): return S.sym_T3p(int(name[4]),int(name[6]))
    if name.startswith('T3_'): return S.sym_T3a(int(name[3]),int(name[5]),int(name[6]))
    raise ValueError(name)
sigma={}
for name,c in atom_terms:
    sym,scalar = atom_sym(name)
    tot=E.fmul(c,scalar)
    for key,cf in sym.items():
        sigma[key]=E.fadd(sigma.get(key,E.ZERO), E.fs(tot,cf))
print("sigma(S)_w3 terms:", len(sigma))
pickle.dump(sigma, open('sigma14.pkl','wb'))
# check: how many distinct tensor keys, and their letter support
from collections import Counter
keyletters = Counter()
for key in sigma:
    keyletters[key]+=1
print("distinct tensor keys:", len(sigma))
# weight-2 part of S (lone L2x,T2) - compute its value for the report
NUM={k:__import__('mpmath').mpf(v) for k,v in pickle.load(open('atoms_num.pkl','rb')).items()}
import mpmath as mp
mp.mp.dps=60
w2val=mp.mpf(0)
for mon,c in S14.items():
    syms=[s for s in dict(mon) if s not in ('l2','la','lu','z2','z3')]
    others=[s for s in dict(mon) if s in ('l2','la','lu','z2','z3')]
    if len(syms)==1 and not others and (syms[0].startswith('L2x') or syms[0].startswith('T2_')):
        w2val += E.fev(c,60)*NUM[syms[0]]
print("weight-2 lone-atom part value:", mp.nstr(w2val,20))
open('progress.log','a').write('PG6 s_B9: sigma(S)_w3 (%d terms) from 44 weight-3 atoms\n'%len(sigma))
