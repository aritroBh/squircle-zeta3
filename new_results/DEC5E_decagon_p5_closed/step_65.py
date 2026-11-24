import sys, pickle
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
QS5 = D.QS5; ZERO, ONE = QS5(0), QS5(1)
P_LET = -D.PHI; Q_LET = D.PHI**2
rewrite = pickle.load(open('../DEC5D_decagon_order5_depth/rewrite_dec_exact.pkl','rb'))
parsed = pickle.load(open('../DEC5D_decagon_order5_depth/parsed_mons_dec.pkl','rb'))
words5 = set()   # words in the |w|+e+f=5 part (Wt5)
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    for n,(p,q) in cu.items():
        for w, gam in rewrite[(n,a,b,c,d)].items():
            if len(w) + e + f == 5 and gam != QS5(0):
                words5.add(w)
print("words in Wt5:", len(words5))
from collections import Counter
def depth(w): return sum(1 for x in w if x != ZERO)
cnt = Counter(depth(w) for w in words5)
print("by depth:", dict(cnt))
# words with letter 1?
w1 = [w for w in words5 if ONE in w]
print("words containing letter 1:", len(w1))
# distinct nonzero-letter multisets
nm = Counter(tuple(sorted(str(x) for x in w if x != ZERO)) for w in words5)
for k, v in sorted(nm.items(), key=lambda kv: -kv[1])[:20]:
    print(v, k)
