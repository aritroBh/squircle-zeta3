
import json
from fractions import Fraction as Fr
D = json.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols.json'))
names = {(1,0,0):'c',(0,1,0):'m',(0,0,1):'p',(-1,0,0):'C',(0,-1,0):'M',(0,0,-1):'P'}
def parse_letter(s):
    digs = []
    i = 0
    while i < len(s):
        if s[i]=='-': digs.append(int(s[i:i+2])); i+=2
        else: digs.append(int(s[i])); i+=1
    return names[tuple(digs)]
def parse_word(w):
    return tuple(parse_letter(p) for p in w.split('|'))
def load(S):
    return {parse_word(w): {int(e): Fr(v) for e,v in cf.items()} for w,cf in S.items()}
R = load(D['R'])
li4 = {k: load(v) for k,v in D['Li4'].items()}
prods = {k: load(v) for k,v in D['prods'].items()}
comps = set()
for S in [R] + list(li4.values()) + list(prods.values()):
    for w, cf in S.items():
        for e in cf:
            comps.add((w, e))
print("total components:", len(comps))
exps = sorted({e for _, e in comps})
print("exponents:", exps)
print("R components:", sum(len(cf) for cf in R.values()))
# R by exponent
for e in exps:
    Re = {w: cf[e] for w, cf in R.items() if e in cf}
    print(f"R exponent {e}: {len(Re)} words:", {''.join(w): str(v) for w,v in sorted(Re.items())})
