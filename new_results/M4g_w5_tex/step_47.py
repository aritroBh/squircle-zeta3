src=open('planner_eval_diag.py').read()
add='''
# dump all needed atom values to atoms_planner.json
need = set()
for nm in d2:
    if nm == 'Lchi8_5': continue
    for f in split(nm): need.add(f)
out = {}
for f in sorted(need):
    v = atom(f)
    out[f] = mp.nstr(v, 460)
out['Lchi8_5'] = mp.nstr(Lchi, 460)
out['N5'] = mp.nstr(N5, 460)
json.dump(out, open('atoms_planner.json','w'))
print('dumped', len(out))
'''
open('planner_eval_diag.py','w').write(src+add)
import subprocess
r=subprocess.run(['nice','-n','10','<project>/.venv/bin/python','planner_eval_diag.py','470'],capture_output=True,text=True,timeout=1200)
print(r.stdout[-300:])
print(r.stderr[-200:])
