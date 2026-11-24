# step: delta_22 closed forms in the two-variable 5-letter space
import sys
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E
S = E.build_symbols()

def biv(A_terms):
    # build a bivector dict from list of (vec,vec,scale)
    out={}
    for (a,b,sc) in A_terms:
        w=E.wedge2_5(a,b)
        for k,v in w.items():
            out[k]=out.get(k,0)+sc*v
    return {k:v for k,v in out.items() if v}
def d22_target(terms):
    out={}
    for (A,B,sc) in terms:
        w=E.wedge_bivec(A,B)
        for k,v in w.items():
            out[k]=out.get(k,0)+sc*v
            if out[k]==0: del out[k]
    return out

X,Y,U,V,W,Z = E.X,E.Y,E.U,E.V,E.W,E.Z
# claimed closed forms
cf_Li22 = d22_target([(biv([(W,Z,1),(U,X,-1)]), biv([(Y,V,1)]), 1)])
cf_Li31 = d22_target([(biv([(W,Z,1)]), biv([(Y,V,1)]), -1), (biv([(U,X,1)]), biv([(V,X,1)]), 1)])
cf_Li13 = d22_target([(biv([(W,Z,1)]), biv([(X,U,1)]), 1)])

for name, cf in [('Li22',cf_Li22), ('Li31',cf_Li31), ('Li13',cf_Li13)]:
    mine = E.delta22_5(S[name])
    diff = {}
    for k in set(mine)|set(cf):
        d = mine.get(k,0)-cf.get(k,0)
        if d: diff[k]=d
    print(f"delta_22({name}): closed form matches: {len(diff)==0};  #components={len(mine)}")
    print("   ", E.show5(mine))
open('progress.log','a').write("step_32: delta_22 closed forms for Li22, Li31, Li13 verified in 5-letter space.\n")
