
base = "<project>/research/08_theory"
import json
d = json.load(open(base + "/G2_M4alpha/hpl_fit_k4_(0,0,0).json"))
print(len(d), "terms")
from collections import Counter
import re
cnt = Counter()
w5 = []
for nm, v in d:
    m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
    if m:
        e = int(m.group(1)); kind = m.group(2)
        w = tuple(int(t) for t in m.group(3).split(',')) if m.group(3) else ()
        cnt[(len(w), kind)] += 1
        if len(w) == 5:
            w5.append((nm, v))
    else:
        cnt[('const', nm)] += 1
print(cnt)
print("\nweight-5 terms:")
for nm, v in w5:
    print(" ", nm, "=", v)
