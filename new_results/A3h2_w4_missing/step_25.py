import sys, json
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup
from evsetup import cc, h_k, sig, a2mp
from flint import acb
import mpmath as mp
mp.mp.dps=130
evsetup.setup(700,1000,800)
BASE='<project>'
rogers=json.load(open(BASE+'/research/08_theory/A3f_level8_relations/weight2_new_relation.json'))['row']

def evalrel_acb(terms,k):
    tot=acb(0)
    for w,c in terms:
        tot += int(c)*h_k(w,k)
    return tot

# 1) shuffle control at weight 2: h(a)h(b)=h(a,b)+h(b,a) for a,b in 1..7 (admissible single letters)
for (a,b) in [(3,5),(1,2),(7,4)]:
    for k in (1,3,5,7):
        d = h_k((a,),k)*h_k((b,),k) - h_k((a,b),k) - h_k((b,a),k)
        print(f'shuffle {a}{b} k={k}: |d|={mp.nstr(abs(a2mp(d)),3)}')
    print()
# 2) five-term with acb summation
print('five-term raw acb:')
for k in (1,3,5,7):
    v=evalrel_acb([(tuple(w),c) for w,c in rogers],k)
    print(f'k={k}: {mp.nstr(abs(a2mp(v)),4)}   (re {mp.nstr(a2mp(v).real,6)})')
