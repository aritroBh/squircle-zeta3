WD='<project>/research/08_theory/FLD1_field_rule_proof'
import pickle
from fractions import Fraction as Fr
import sys
sys.path.insert(0,WD)
from hpl import coeffs_exact
store = pickle.load(open(WD+'/hpl_constructed_K5.pkl','rb'))
def series_mine(terms,lo,hi):
    out={}
    for (e,w,q) in terms:
        if len(w)==0:
            if lo<=e<=hi: out[e]=out.get(e,Fr(0))+q
            continue
        ce=coeffs_exact(w,hi-e+2)
        for j in range(lo,hi+1):
            idx=j-e
            if 1<=idx<len(ce): out[j]=out.get(j,Fr(0))+q*ce[idx]
    return out
def Fclosed(ca, ct, hi):
    out={}
    for e,q in ca.items():
        for j in range(1,hi-e+2,2):  # artanh c^j/j, odd j, j+e<=hi
            out[j+e]=out.get(j+e,Fr(0))+q*Fr(1,j)
    for e,q in ct.items(): out[e]=out.get(e,Fr(0))+q
    return out
hi=40
F1 = Fclosed({1:Fr(7,16),-1:Fr(-5,8),-3:Fr(3,16)}, {0:Fr(9,16),-2:Fr(-3,16)}, hi)
F2 = Fclosed({1:Fr(-1,16),-1:Fr(-1,8),-3:Fr(3,16)}, {0:Fr(1,16),-2:Fr(-3,16)}, hi)
m4=series_mine(store[(4,(0,0,1,0))],-4,hi); m22=series_mine(store[(4,(2,0,0,0))],-4,hi)
for name,g,mine in [('z4 F1',F1,m4),('z2^2 F2',F2,m22)]:
    diff={j:g.get(j,Fr(0))-mine.get(j,Fr(0)) for j in range(2,hi+1)}
    diff={j:v for j,v in diff.items() if v!=0}
    print(f'{name}: agree (c^2..c^{hi}): {len(diff)==0}', ('' if not diff else list(diff.items())[:5]))