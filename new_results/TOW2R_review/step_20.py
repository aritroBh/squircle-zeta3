import json, os, re
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
print('=== planner/prove_adj34.out ===')
print(open(os.path.join(tow2,'planner/prove_adj34.out')).read())
print('=== planner/prove_tops.out ===')
print(open(os.path.join(tow2,'planner/prove_tops.out')).read())
tr = json.load(open(os.path.join(tow2, 'transcript.json')))
txt = json.dumps(tr)
for pat in [r'adjacent \(\d,\d\): remainder\*DEN == 0 : \w+', r'E3 blindness delta=\d: numerator == 0 : \w+',
            r'S\^\d/S identity[^\\]*', r'l=\d: total built[^\\]*']:
    print(pat, '->', re.findall(pat, txt)[:10])
