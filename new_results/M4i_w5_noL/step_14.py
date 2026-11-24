import json
base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
r5 = json.loads(open(m4f + '/R5_best.json').read())
print('R5_best keys:', list(r5.keys()))
terms = r5['terms']
print('n terms:', len(terms))
print('Lchi8_5:', terms.get('Lchi8_5'))
print('z5:', terms.get('z5'))
li5 = [k for k in terms if k.startswith('Li5')]
print('Li5 terms (%d):' % len(li5))
for k in li5: print('  ', k, terms[k])
print('sample products:', [k for k in terms if '*' in k][:6])
print()
print('Lchi8_5.json:', open(m4f + '/Lchi8_5.json').read()[:300])
print()
print('W5_5400.json:', open(m4f + '/W5_5400.json').read()[:300])
print()
print('progress.log:'); print(open(m4f + '/progress.log').read())
import os
print('planner dir:', os.listdir(m4f + '/planner'))
