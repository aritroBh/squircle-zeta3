
import json
from fractions import Fraction as Fr
D = json.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols.json'))
def lw(s):
    # letter word string "010|..." -> readable
    names = {(1,0,0):'c',(0,1,0):'m',(0,0,1):'p',(-1,0,0):'C',(0,-1,0):'M',(0,0,-1):'P'}
    out = []
    i = 0
    digs = []
    while i < len(s):
        if s[i]=='-': digs.append(int(s[i:i+2])); i+=2
        else: digs.append(int(s[i])); i+=1
    return names[tuple(digs)]
def parse_word(w):
    parts = w.split('|')
    return tuple(lw(p) for p in parts)
R = {parse_word(w): {int(e): Fr(v) for e,v in cf.items()} for w,cf in D['R'].items()}
DT = {parse_word(w): {int(e): Fr(v) for e,v in cf.items()} for w,cf in D['DT'].items()}
print("DT terms:")
for w,cf in sorted(DT.items()):
    print("  ", ''.join(w), {e: str(v) for e,v in sorted(cf.items())})
print()
print("R = S(Omega-DT): exponents present:", sorted({e for cf in R.values() for e in cf}))
for w,cf in sorted(R.items()):
    print("  ", ''.join(w), {e: str(v) for e,v in sorted(cf.items())})
