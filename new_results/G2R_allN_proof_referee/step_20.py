import json, re
G2 = "<project>/research/08_theory/G2_M4alpha/"
d = json.load(open(G2+"hpl_fit_k4_(0,0,0).json"))
from collections import Counter
cnt = Counter()
for nm, v in d:
    m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
    w = tuple(int(t) for t in m.group(3).split(','))
    cnt[len(w)] += 1
print(cnt, "total", sum(cnt.values()))
print("weight-2 terms:", [(nm,v) for nm,v in d if nm.count(',')==0 and 'He' in nm])
