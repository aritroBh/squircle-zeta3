import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
cands = json.load(open(c+'/out/d13/candidates_cap.json'))
print('keys:', list(cands.keys()))
v = cands[list(cands.keys())[0]]
print(json.dumps({k: cands[k] for k in list(cands)[:1]}, indent=1)[:500])
# maybe there's a 'list' key
for k in cands:
    print(k, type(cands[k]), (len(cands[k]) if hasattr(cands[k],'__len__') else ''))
