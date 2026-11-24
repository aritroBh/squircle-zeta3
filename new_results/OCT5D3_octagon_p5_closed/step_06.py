import pickle
base = "<project>/research/08_theory/"

# parsed_mons
parsed = pickle.load(open(base+"OCT5D_octagon_order5_depth/parsed_mons.pkl","rb"))
print(type(parsed), len(parsed))
print("first entry:", parsed[0])
# count weight-5 monomials
w5 = [m for m in parsed if sum(m[0])==5]
print("weight-5 monomials:", len(w5))
from collections import Counter
degcount = Counter(sum(m[0]) for m in parsed)
print("degree counts:", dict(sorted(degcount.items())))

# p5_800
p5 = pickle.load(open(base+"OCT5D_octagon_order5_depth/p5_800.pkl","rb"))
print("\np5_800 type:", type(p5))
if isinstance(p5, dict):
    print(p5.keys())
    for k,v in p5.items():
        print(k, type(v), str(v)[:100])
else:
    print(str(p5)[:500])
