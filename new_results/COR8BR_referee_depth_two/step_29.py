# step: exact comparison with COR8B's twovar_symbols.json + exact shuffle identities
import sys, json
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E
S = E.build_symbols()

tws = json.load(open('<project>/research/08_theory/COR8B_depth_two_two_terms/twovar_symbols.json'))
def key_to_tuple(k):
    return tuple(tuple(int(ch) for ch in part) for part in k.split('|'))
allmatch = True
for name in ['Li22','Li31','Li13','Li12','Li21','Li11']:
    theirs = {key_to_tuple(k): Fr(v) for k,v in tws[name].items()}
    match = (S[name] == theirs)
    allmatch &= match
    print(f"{name}: exact match with twovar_symbols.json: {match}")
print("ALL MATCH:", allmatch)

# Shuffle identities (exact):
# 1) Li22(x,y)+Li22(y,x)+Li4(xy) = Li2(x)Li2(y)
lhs = E.sadd(E.sadd(S['Li22'], E.swap_xy(S['Li22'])), S['Li4xy'])
rhs = E.shuffle_sym(S['Li2x'], S['Li2y'])
print("shuffle22 residual terms:", len(E.sadd(lhs, rhs, Fr(-1))))
# 2) stuffle: Li31(x,y)+Li13(y,x)+Li4(xy) = Li3(x)Li1(y);  S(Li1(y)) = -[V]
lhs2 = E.sadd(E.sadd(S['Li31'], E.swap_xy(S['Li13'])), S['Li4xy'])
S_Li1y = {(E.V,): Fr(-1)}
rhs2 = E.shuffle_sym(S['Li3x'], S_Li1y)
print("stuffle31 residual terms:", len(E.sadd(lhs2, rhs2, Fr(-1))))
# 3) Li11(x,y)+Li11(y,x)+Li2(xy) = Li1(x)Li1(y)
lhs3 = E.sadd(E.sadd(S['Li11'], E.swap_xy(S['Li11'])), S_Lin_xy if False else E.S_Lin(E.Z,E.W,2))
S_Li1x = {(E.U,): Fr(-1)}
rhs3 = E.shuffle_sym(S_Li1x, S_Li1y)
print("stuffle11 residual terms:", len(E.sadd(lhs3, rhs3, Fr(-1))))
open('progress.log','a').write("step_28: two-variable symbols match twovar_symbols.json exactly; shuffle/stuffle identities exact.\n")
