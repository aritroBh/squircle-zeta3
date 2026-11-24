import os, json
dst = '<project>/research/08_theory/M4i_w5_noL'
for f in ['step_02.py','step_03.py','step_04.py','step_05.py']:
    p = os.path.join(dst,f)
    print('='*20, f, os.path.getsize(p))
    print(open(p).read()[:800])
tr = json.loads(open(dst+'/transcript.json').read())
print('transcript type:', type(tr), 'len:', len(tr) if hasattr(tr,'__len__') else '')
