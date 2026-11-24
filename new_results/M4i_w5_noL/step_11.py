import json
base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
oc = json.loads(open(m4f + '/overcomplete_1200.json').read())
print('keys:', list(oc.keys()))
print('kept:', oc['kept'])
print('extras:', oc['extras'])
print('n deps:', len(oc['deps']))
print('dep[0] keys:', list(oc['deps'][0].keys()))
print('dep[0] kept len:', len(oc['deps'][0]['kept']))
print('dep[0] vec:', oc['deps'][0]['vec'])
print('dep[0] resid:', oc['deps'][0]['resid'])
print()
tr = oc['target_rels'][0]
print('target rel maxcoef:', tr['maxcoef'], 'resid:', tr['resid'])
print('target rel vec:', tr['vec'])
