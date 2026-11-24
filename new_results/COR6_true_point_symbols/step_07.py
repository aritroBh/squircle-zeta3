base = "<project>/research/"
import json
D = json.load(open(base + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
print(type(D), len(D))
for k, v in D[:15]:
    print(k, "->", v)
print("...")
# weights present
from collections import Counter
cnt = Counter()
for k, v in D:
    e, w = k.split(' ', 1)
    w = eval(w[1:])
    cnt[len(w)] += 1
print(cnt)
