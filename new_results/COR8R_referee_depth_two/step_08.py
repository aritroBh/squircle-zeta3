import json, os, re
base = '<project>'
data = json.load(open(os.path.join(base,'research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json')))
from collections import Counter
# parse terms: "c^E H[w1, w2, ...]" with coefficient
terms = []
for name, coeff in data:
    m = re.match(r'c\^(-?\d+)\s*H\[([^\]]*)\]', name)
    E = int(m.group(1))
    w = tuple(int(x) for x in m.group(2).split(','))
    terms.append((E, w, int(coeff)))
cnt = Counter(len(w) for _,w,_ in terms)
print("weights present:", dict(cnt))
print("exponents E:", sorted(set(E for E,_,_ in terms)))
wt4 = [(E,w,c) for E,w,c in terms if len(w)==4]
print("num weight-4 terms:", len(wt4))
for t in wt4[:15]: print(t)
