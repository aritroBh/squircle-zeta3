
import subprocess, os
REV = 'research/08_theory/D13r_lemma_review'
PY='<home>)/.venv/bin/python'
for f in ['step_31.py']:
    r = subprocess.run([PY, os.path.join(REV,f)], capture_output=True, text=True, cwd=REV, timeout=600)
    print("="*20, f, "="*20)
    print(r.stdout)
    print(r.stderr[-2000:] if r.stderr else "")
