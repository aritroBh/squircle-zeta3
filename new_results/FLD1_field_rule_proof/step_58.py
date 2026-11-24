WD='<project>/research/08_theory/FLD1_field_rule_proof'
import json, pickle, re
from fractions import Fraction as Fr
import sys
sys.path.insert(0,WD)
from hpl import coeffs_exact
G1='<project>/research/10_frontier/G1_zeta3_any_angle/'
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
def load_g1(path):
    out=[]
    for nm,v in json.load(open(path)):
        v=Fr(v)
        if nm=='1': out.append((0,(),v)); continue
        if nm=='c': out.append((1,(),v)); continue
        mm=re.match(r'c\^(-?\d+) H\[(.*)\]', nm); out.append((int(mm.group(1)),tuple(int(t) for t in mm.group(2).split(',')),v))
    return out
hi=40
# z4 F1, z2^2 F2 (restrict to c^2..c^40)
def artanh_series(hi): return {j: Fr(1,j) for j in range(1,hi+1,2)}
def Fclosed(ca, ct, hi):
    at=artanh_series(hi+2); out={}
    for e,q in ca.items():
        for j,v in at.items():
            if j+e<=hi: out[j+e]=out.get(j+e,Fr(0))+q*v
    for e,q in ct.items(): out[e]=out.get(e,Fr(0))+q
    return out
F1 = Fclosed({1:Fr(7,16),-1:Fr(-5,8),-3:Fr(3,16)}, {0:Fr(9,16),-2:Fr(-3,16)}, hi)
F2 = Fclosed({1:Fr(-1,16),-1:Fr(-1,8),-3:Fr(3,16)}, {0:Fr(1,16),-2:Fr(-3,16)}, hi)
m4=series_mine(store[(4,(0,0,1,0))],-4,hi); m22=series_mine(store[(4,(2,0,0,0))],-4,hi)
for name,g,mine in [('z4 F1',F1,m4),('z2^2 F2',F2,m22)]:
    diff={j:g.get(j,Fr(0))-mine.get(j,Fr(0)) for j in range(2,hi+1)}
    diff={j:v for j,v in diff.items() if v!=0}
    print(f'{name}: agree (c^2..c^{hi}): {len(diff)==0}', ('' if not diff else list(diff.items())[:5]))
# z3: check my z3 == -(G1 B3 z2)
g1z2 = load_g1(G1+'hpl_fit_k3_(1,0,0).json')
gz = series_mine(g1z2,-4,hi)
mine_z3 = series_mine(store[(4,(0,1,0,0))],-4,hi)
for sgn,lab in [(1,'+B3z2'),(-1,'-B3z2')]:
    diff={j:sgn*gz.get(j,Fr(0))-mine_z3.get(j,Fr(0)) for j in range(-4,hi+1)}
    diff={j:v for j,v in diff.items() if v!=0}
    print(f'my k4 z3 == {lab}: {len(diff)==0}', ('' if not diff else list(diff.items())[:5]))