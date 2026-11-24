
import json
from fractions import Fraction as Fr
D = json.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols.json'))
def parse_letter(s):
    digs = []
    i = 0
    while i < len(s):
        if s[i]=='-': digs.append(int(s[i:i+2])); i+=2
        else: digs.append(int(s[i])); i+=1
    return tuple(Fr(x) for x in digs)
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
comps = sorted(comps, key=str)
print("total components:", len(comps), " unknowns:", len(li4) + len(prods))
exps = sorted({e for _, e in comps})
print("exponents:", exps)
for e in exps:
    Re = {w: cf[e] for w, cf in R.items() if e in cf}
    def rs(l):
        b = {(1,0,0):'c',(0,1,0):'m',(0,0,1):'p'}
        parts=[]
        for j,nm in enumerate('cmp'):
            if l[j]: parts.append(nm + (f'^{l[j]}' if l[j]!=1 else ''))
        return '*'.join(parts)
    print(f"R e={e}:", { '*'.join(rs(l) for l in w): str(v) for w,v in sorted(Re.items(), key=str)})
import pickle
pickle.dump({'R':R,'li4':li4,'prods':prods,'comps':comps}, open('solve_data.pkl','wb'))
