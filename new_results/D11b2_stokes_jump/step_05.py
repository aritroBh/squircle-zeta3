import json
pj = json.load(open('pade1_9090.json'))
print('pade1_9090 keys:', list(pj.keys()), 'len num:', len(pj['num']), 'len den:', len(pj['den']))
print('part1_final:', json.load(open('part1_final.json')))
d = json.load(open('step2_num450_N180.json'))
print('step2 keys:', list(d.keys()), 'N:', len(d['f450']))
print(open('../D11b_stokes_numeric/result.md').read()[:3000])
