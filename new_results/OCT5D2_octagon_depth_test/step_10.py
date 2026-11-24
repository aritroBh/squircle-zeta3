import pickle
base = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(base+"parsed_mons.pkl","rb"))
print(type(parsed), len(parsed))
print("first 3 entries:")
for mon in parsed[:3]:
    print(mon)
print("...")
# structure: list of (mon_tuple, coeff_dict)
mons = [m for m,c in parsed]
from collections import Counter
print("degree counts:", Counter(sum(m) for m in mons))
# class counts at weight 5
def cls(m):
    a,b,c,d,e,f = m
    if c>=1 and d>=1: return 'mixed'
    if c>=1 or d>=1: return 'single'
    return 'local'
print("weight-5 class counts:", Counter(cls(m) for m in mons if sum(m)==5))
# coefficient dict structure
print("sample coeff dicts:")
for m,c in parsed[:5]:
    print(m, c)
