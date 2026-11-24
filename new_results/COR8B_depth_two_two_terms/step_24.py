
import json
d = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
import re
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0])
    e = int(m.group(1)); rest = m.group(2).strip()
    if rest.startswith('H'):
        w = [int(x) for x in rest[2:-1].split(',')]
        return e, w
    else:
        return e, rest
from collections import Counter
cnt = Counter()
for t in d:
    e, w = parse(t)
    cnt[(len(w) if isinstance(w,list) else 'const')] += 1
print(cnt)
# print all weight-4 terms
print("\nweight-4 terms:")
for t in d:
    e, w = parse(t)
    if isinstance(w, list) and len(w) == 4:
        print(f"c^{e} H{w}  coeff {t[1]}")
print("\nnon-weight-4 terms:")
for t in d:
    e, w = parse(t)
    if not (isinstance(w, list) and len(w) == 4):
        print(t[0], "coeff", t[1])
