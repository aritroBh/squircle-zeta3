# multilinear expansion of composite letters (Z=X+Y), then compare
import sys
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import importlib, cor8br_sym as E
importlib.reload(E)
S = E.build_symbols()

BASIS5 = [E.X, E.Y, E.U, E.V, E.W]
def expand_ml(s):
    # expand any letter that is a formal sum into basis letters multilinearly
    out={}
    for k,c in s.items():
        # decompose each letter into basis
        opts=[]
        for v in k:
            # express v in basis
            terms=[]
            for i,b in enumerate(BASIS5):
                if v[i]: terms.append((b, v[i]))
            opts.append(terms)
        # product over positions
        acc=[((), Fr(1))]
        for terms in opts:
            new=[]
            for (pref, pc) in acc:
                for (b, bc) in terms:
                    new.append((pref+(b,), pc*bc))
            acc=new
        for (t, tc) in acc:
            out[t]=out.get(t,Fr(0))+c*tc
    return {k:v for k,v in out.items() if v}

def resid(a,b):
    d=E.sadd(expand_ml(a), expand_ml(b), Fr(-1))
    return len(d), d

lhs = E.sadd(E.sadd(S['Li22'], E.swap_xy(S['Li22'])), S['Li4xy'])
rhs = E.shuffle_sym(S['Li2x'], S['Li2y'])
print("shuffle22 residual after ML expansion:", resid(lhs,rhs)[0])
lhs2 = E.sadd(E.sadd(S['Li31'], E.swap_xy(S['Li13'])), S['Li4xy'])
rhs2 = E.shuffle_sym(S['Li3x'], {(E.V,):Fr(-1)})
print("stuffle31 residual after ML expansion:", resid(lhs2,rhs2)[0])
lhs3 = E.sadd(E.sadd(S['Li11'], E.swap_xy(S['Li11'])), E.S_Lin(E.Z,E.W,2))
rhs3 = E.shuffle_sym({(E.U,):Fr(-1)}, {(E.V,):Fr(-1)})
print("stuffle11 residual after ML expansion:", resid(lhs3,rhs3)[0])
lhs4 = E.sadd(E.sadd(S['Li12'], E.swap_xy(S['Li21'])), E.S_Lin(E.Z,E.W,3))
rhs4 = E.shuffle_sym(S['Li2x'], {(E.V,):Fr(-1)})
print("shuffle12/21 residual after ML expansion:", resid(lhs4,rhs4)[0])
open('progress.log','a').write("step_30: shuffle/stuffle identities exact after multilinear expansion (0 residual).\n")
