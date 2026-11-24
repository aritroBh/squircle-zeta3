# step: build two-variable symbols from my recursion; compare with COR8B twovar_symbols.json; exact shuffle checks
import sys, json
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E

S = E.build_symbols()
for name in ['Li11','Li12','Li21','Li22','Li31','Li13']:
    print(name, "terms:", len(S[name]))

# compare with twovar_symbols.json (basis order there: 0=X,1=Y,2=U,3=V,4=W -- same as mine)
tws = json.load(open('<project>/research/08_theory/COR8B_depth_two_two_terms/twovar_symbols.json'))
def key_to_tuple(k):
    vecs = tuple(tuple(int(ch) for ch in part) for part in k.split('|'))
    return vecs
for name in ['Li22','Li31','Li13','Li12','Li21','Li11']:
    theirs = {key_to_tuple(k): Fr(v) for k,v in tws[name].items()}
    mine = S[name]
    match = (mine == theirs)
    print(f"{name}: mine {len(mine)} terms, theirs {len(theirs)} terms, exact match: {match}")
    if not match:
        only_mine = {k:v for k,v in mine.items() if theirs.get(k)!=v}
        only_theirs = {k:v for k,v in theirs.items() if mine.get(k)!=v}
        print("  diff mine:", only_mine)
        print("  diff theirs:", only_theirs)
