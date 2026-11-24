
import pickle
base = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
ITAB = pickle.load(open(base+"itab110.pkl","rb"))
ks = list(ITAB.keys())
print(len(ks), ks[:5])
parsed = pickle.load(open(base+"parsed_mons.pkl","rb"))
print(len(parsed), parsed[0])
# weight-5 monomials
w5 = [(mon,cu) for mon,cu in parsed if sum(mon)==5]
print("weight-5 monomials:", len(w5))
import collections
cnt = collections.Counter()
for mon,cu in w5:
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    cnt[cls]+=1
print(cnt)
