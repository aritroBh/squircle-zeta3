import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
cc = json.load(open(c+'/out/d13/candidates_cap.json'))
print('meta:', {k:cc[k] for k in ['D','B','split','inert','names','S','cap','kernel_dim']})
c0 = cc['cands'][0]
print(json.dumps(c0, indent=1))
print('...')
print(json.dumps(cc['cands'][274], indent=1)[:400])
