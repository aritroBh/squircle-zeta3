import sys
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E
S = E.build_symbols()
BASIS5 = [E.X, E.Y, E.U, E.V, E.W]
def expand_ml(s):
    out={}
    for k,c in s.items():
        opts=[]
        for v in k:
            terms=[(b,v[i]) for i,b in enumerate(BASIS5) if v[i]]
            opts.append(terms)
        acc=[((), Fr(1))]
        for terms in opts:
            acc=[(pref+(b,), pc*bc) for (pref,pc) in acc for (b,bc) in terms]
        for (t,tc) in acc:
            out[t]=out.get(t,Fr(0))+c*tc
    return {k:v for k,v in out.items() if v}
# correct identity: Li21(x,y) + Li12(y,x) + Li3(xy) = Li2(x)Li1(y)
lhs = E.sadd(E.sadd(S['Li21'], E.swap_xy(S['Li12'])), E.S_Lin(E.Z,E.W,3))
rhs = E.shuffle_sym(S['Li2x'], {(E.V,):Fr(-1)})
d = E.sadd(expand_ml(lhs), expand_ml(rhs), Fr(-1))
print("shuffle Li2(x)Li1(y) residual:", len(d))
open('progress.log','a').write("step_31: correct weight-3 shuffle Li21(x,y)+Li12(y,x)+Li3(xy)=Li2(x)Li1(y) exact (0 residual); earlier pairing was my own guess, not COR8B's claim.\n")
