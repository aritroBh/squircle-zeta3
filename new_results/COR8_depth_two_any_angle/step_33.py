
import json
from fractions import Fraction as Fr
D = json.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols.json'))
def parse_letter(s):
    digs = []
    i = 0
    while i < len(s):
        if s[i]=='-': digs.append(int(s[i:i+2])); i+=2
        else: digs.append(int(s[i])); i+=1
    return tuple(digs)
names = {(1,0,0):'c',(0,1,0):'m',(0,0,1):'p',(-1,0,0):'C',(0,-1,0):'M',(0,0,-1):'P',(0,0,0):'0'}
def parse_word(w):
    return tuple(names[parse_letter(p)] for p in w.split('|'))
DT = {parse_word(w): {int(e): Fr(v) for e,v in cf.items()} for w,cf in D['DT'].items()}
print("DT terms (raw, with zero letters):")
for w,cf in sorted(DT.items()):
    print("  ", ''.join(w), {e: str(v) for e,v in sorted(cf.items())})
