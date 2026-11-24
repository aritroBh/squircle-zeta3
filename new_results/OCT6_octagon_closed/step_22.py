import pickle, sympy as sp, mpmath as mp
import mixed_engine as ME
mons = pickle.load(open('mons_simple.pkl','rb'))

c1log = -ME.l2/2 - ME.lu
c2log = -ME.lu

S_local = sp.Integer(0); S_s1 = sp.Integer(0); S_s2 = sp.Integer(0); S_mix = sp.Integer(0)
for (a,b,c1,c2,e1,e2), terms in mons:
    clog = c1log**e1 * c2log**e2
    for n,q in terms:
        J = ME.INT(n,a,b,c1,c2)
        piece = sp.expand(clog * q * J)
        if c1==0 and c2==0: S_local += piece
        elif c2==0: S_s1 += piece
        elif c1==0: S_s2 += piece
        else: S_mix += piece
S_tot = sp.expand(S_local + S_s1 + S_s2 + S_mix)
S_mix = sp.expand(S_mix)
print("assembly done; S_mix terms:", len(sp.Add.make_args(S_mix)), " S_tot terms:", len(sp.Add.make_args(S_tot)))
pickle.dump({'S_local':sp.expand(S_local),'S_s1':sp.expand(S_s1),'S_s2':sp.expand(S_s2),'S_mix':S_mix,'S_tot':S_tot}, open('S_assembly.pkl','wb'))
print("saved S_assembly.pkl")
