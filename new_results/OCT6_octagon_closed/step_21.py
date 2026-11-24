import pickle, sympy as sp, mpmath as mp
import mixed_engine as ME
from mixed_engine import INT, NAMED, l2, lu, c1 if False else None
mons = pickle.load(open('mons_simple.pkl','rb'))
r2 = sp.sqrt(2); t = sp.symbols('t')

# constant logs c1log = -l2/2 - lu, c2log = -lu
c1log = -ME.l2/2 - ME.lu
c2log = -ME.lu

# assemble S (so that int E3 = u*S); separate mixed part for its own check
S_local = sp.Integer(0); S_s1 = sp.Integer(0); S_s2 = sp.Integer(0); S_mix = sp.Integer(0)
for (a,b,c1,c2,e1,e2), terms in mons:
    clog = c1log**e1 * c2log**e2
    for n,q in terms:
        J = INT(n,a,b,c1,c2)
        piece = sp.expand(clog * q * J)
        if c1==0 and c2==0: S_local += piece
        elif c2==0: S_s1 += piece
        elif c1==0: S_s2 += piece
        else: S_mix += piece
S_tot = sp.expand(S_local + S_s1 + S_s2 + S_mix)
S_mix = sp.expand(S_mix)
print("assembly done; S_mix terms:", len(sp.Add.make_args(S_mix)), " S_tot terms:", len(sp.Add.make_args(S_tot)))
pickle.dump({'S_local':S_local,'S_s1':S_s1,'S_s2':S_s2,'S_mix':S_mix,'S_tot':S_tot}, open('S_assembly.pkl','wb'))
