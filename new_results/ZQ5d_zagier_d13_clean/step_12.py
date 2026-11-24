import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
k5 = json.load(open(c+'/out/d13/kernel5_cap.json'))
print('names:', k5['names'], 'd:', k5['d'], 'rank:', k5['rank'])
K = k5['K']
print('K rows:', len(K), 'len row0:', len(K[0]))
print('K[0][:20]:', K[0][:20])
comps = k5['comps']
print('comps[0]:', comps[0])
cands = json.load(open(c+'/out/d13/candidates_cap.json'))
print('cands:', type(cands), len(cands))
print(json.dumps(cands[0], indent=1)[:800])
