import subprocess, os
r=subprocess.run(['nice','-n','10','/opt/homebrew/bin/tectonic','-X','compile','fragment_test.tex'],capture_output=True,text=True,timeout=600)
w=[l for l in (r.stdout+r.stderr).splitlines() if 'warning' in l.lower() or 'error' in l.lower()]
print('fragment_test:', r.stdout.strip().splitlines()[-2].strip(), '| warnings:', w if w else 'NONE')
print('--- paper_fragment.tex (current, verbatim) ---')
print(open('paper_fragment.tex').read())
print('--- fragment_test.tex ---')
print(open('fragment_test.tex').read())
for f in sorted(os.listdir('.')):
    if not f.startswith('transcript') and not f.startswith('step_') and f!='bisect2.py' and f!='planner_eval_diag.py' and f!='atoms_planner.json':
        print(f, os.path.getsize(f))
