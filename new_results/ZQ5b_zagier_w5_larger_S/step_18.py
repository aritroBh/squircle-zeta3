import os, json
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
# Check candidate counts and h distribution for D=12, 21, 17
for D in (12, 21, 17, 8, 5):
    J = json.load(open(f'../ZQ5_zagier_weight5/out/d{D}/candidates.json'))
    hs = [c['h'] for c in J['cands']]
    import collections
    print(D, len(hs), 'h dist:', dict(sorted(collections.Counter(hs).items())))
