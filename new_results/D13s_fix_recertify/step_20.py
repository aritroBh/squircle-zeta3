
import subprocess, os
REV = 'research/08_theory/D13r_lemma_review'
PY='<home>)/.venv/bin/python'
r = subprocess.run([PY, os.path.join(REV,'step_32.py')], capture_output=True, text=True, cwd=REV, timeout=600)
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else "")
