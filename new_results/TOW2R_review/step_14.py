import os, re
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
# grep for value identity proofs: look for "V(3)", "V(4)", "value", "(3,0)", "(4,0)" in step files
hits = []
for fn in sorted(os.listdir(tow2)):
    if not fn.endswith('.py'): continue
    src = open(os.path.join(tow2, fn)).read()
    for kw in ['V(3)', 'V(4)', 'value identity', '(3,0)', '(4,0)', 'V3', 'V4', 'slope', 'value']:
        if kw in src:
            hits.append((fn, kw))
from collections import defaultdict
byfile = defaultdict(list)
for fn, kw in hits: byfile[fn].append(kw)
for fn in sorted(byfile): print(fn, sorted(set(byfile[fn])))
