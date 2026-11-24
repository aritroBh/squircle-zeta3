WD='<project>/research/08_theory/FLD1_field_rule_proof'
# quick: compare my k4 z4, z2^2 parts to G2's closed forms F1,F2 (b4_pi4_part.json) and z3 to G1's B3 z2 part
import json, pickle, re
from fractions import Fraction as Fr
import sys
sys.path.insert(0,WD)
from hpl import coeffs_exact
G2='<project>/research/08_theory/G2_M4alpha/'
G1='<project>/research/10_frontier/G1_zeta3_any_angle/'
print("b4_pi4_part.json:", json.load(open(G2+'b4_pi4_part.json')))
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
# G2 F1, F2 as Laurent series in c (artanh c = sum c^{2j+1}/(2j+1))
def artanh_series(hi):
    return {j: Fr(1,j) for j in range(1,hi+1,2)}  # c^j/j for odd j
def Fclosed(coef_artanh, const_terms, hi):
    # coef_artanh: dict e->q meaning (sum_e q c^e)*artanh(c); const_terms: dict e->q bare
    at = artanh_series(hi+4)
    out={}
    for e,q in coef_artanh.items():
        for j,v in at.items():
            out[j+e]=out.get(j+e,Fr(0))+q*v
    for e,q in const_terms.items():
        out[e]=out.get(e,Fr(0))+q
    return out
hi=40
# F1 = (7c/16 - 5/(8c) + 3/(16c^3)) artanh + 9/16 - 3/(16c^2)
F1 = Fclosed({1:Fr(7,16),-1:Fr(-5,8),-3:Fr(3,16)}, {0:Fr(9,16),-2:Fr(-3,16)}, hi)
# F2 = (-c/16 - 1/(8c) + 3/(16c^3)) artanh + 1/16 - 3/(16c^2)
F2 = Fclosed({1:Fr(-1,16),-1:Fr(-1,8),-3:Fr(3,16)}, {0:Fr(1,16),-2:Fr(-3,16)}, hi)
m4 = series_mine(store[(4,(0,0,1,0))],-4,hi)  # z4
m22 = series_mine(store[(4,(2,0,0,0))],-4,hi) # z2^2
for name,g,mine in [('z4 F1',F1,m4),('z2^2 F2',F2,m22)]:
    diff={j:g.get(j,Fr(0))-mine.get(j,Fr(0)) for j in set(g)|set(mine) if 0<=j}
    diff={j:v for j,v in diff.items() if v!=0 and j>=2}
    print(f'{name}: my k4 part vs G2 closed form (even coeffs c^2..c^{hi}) agree: {len(diff)==0}', ('' if not diff else list(diff.items())[:5]))
# z3 part vs G1 B3 z2 part (hpl_fit_k3_(1,0,0).json), using He convention? it's G1 json (H not He)
def load_g1(path):
    out=[]
    for nm,v in json.load(open(path)):
        v=Fr(v)
        if nm=='1': out.append((0,(),v)); continue
        if nm=='c': out.append((1,(),v)); continue
        mm=re.match(r'c\^(-?\d+) H\[(.*)\]', nm); out.append((int(mm.group(1)),tuple(int(t) for t in mm.group(2).split(',')),v))
    return out
g1z2 = load_g1(G1+'hpl_fit_k3_(1,0,0).json')
def series_g1(terms,lo,hi):
    return series_mine(terms,lo,hi)
mine_z3 = series_mine(store[(4,(0,1,0,0))],-4,hi)
g_z3 = series_g1(g1z2,-4,hi)
diff={j:g_z3.get(j,Fr(0))-mine_z3.get(j,Fr(0)) for j in set(g_z3)|set(mine_z3)}
diff={j:v for j,v in diff.items() if v!=0}
print('z3 part: mine vs G1 B3-z2 (G2 says k4 z3 == B3 z2): agree:', len(diff)==0, ('' if not diff else list(diff.items())[:5]))