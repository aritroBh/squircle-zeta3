import json, subprocess
d=json.load(open('../M4f_w5_extras/R5_best.json'))
print('keys:', list(d.keys())[:8])
t=d.get('terms', d)
print('n terms:', len(t))
print('sample:', list(t.items())[:3])
li5=[k for k in t if k.startswith('Li5(')]
print('Li5 args:', sorted(li5))
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'gen_supplement.py','../M4f_w5_extras/R5_best.json','supplement_W5_best.tex'],capture_output=True,text=True)
print(r.stdout.strip().splitlines()[-3:]); print(r.stderr[-400:])
