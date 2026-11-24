
import subprocess, os
REV = 'research/08_theory/D13r_lemma_review'
PY='<home>)/.venv/bin/python'
for f in ['step_29.py','step_30.py']:
    r = subprocess.run([PY, os.path.join(REV,f)], capture_output=True, text=True, cwd=REV, timeout=900)
    print("="*20, f, "="*20); print(r.stdout); print(r.stderr[-800:] if r.stderr else "")
