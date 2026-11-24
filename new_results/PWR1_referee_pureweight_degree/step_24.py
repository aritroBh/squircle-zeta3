import json
js = json.load(open('../D16_pure_weight/tau_series_2d.json'))
for k in ['P_over_sqrt2','logP','B2']:
    v = js[k]
    print(k, type(v))
    s = str(v)
    print(s[:800])
    print('---')
