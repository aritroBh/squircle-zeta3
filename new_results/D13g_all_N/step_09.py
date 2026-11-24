
import json, os
BASE='research/08_theory'
# truth data
p = os.path.join(BASE,'D13f_sweep','mk_deg9.tsv')
print(open(p).read()[:800])
print('...')
# step_34 sharp pieces
q = os.path.join(BASE,'D13f_sweep','cert_single_N.json')
d = json.load(open(q))
for N,v in d.items():
    print(N, {k: v[k][:24] for k in ['C','W1','W2','W3','TM','TR','TLm','TLf','TL0','D','rho']})
