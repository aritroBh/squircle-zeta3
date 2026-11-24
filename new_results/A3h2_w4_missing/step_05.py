import json
base='<project>'
d=json.load(open(base+'/research/08_theory/A3h_weight4_target/T3_quotient_nf.json'))
print('p =',d['p'])
print('freewords (n=%d):'%len(d['freewords']))
print(d['freewords'])
print()
print('P4_nf (first 30):', d['P4_nf'][:30] if isinstance(d['P4_nf'],list) else type(d['P4_nf']))
print('len Drat_nf', len(d['Drat_nf']), 'len Dsq_nf', len(d['Dsq_nf']), 'len P4_nf', len(d['P4_nf']))
print('Drat_nf[:8]', d['Drat_nf'][:8])
