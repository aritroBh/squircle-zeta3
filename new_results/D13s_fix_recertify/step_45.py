
import subprocess, os
REV = 'research/08_theory/D13r_lemma_review'
PY='<home>)/.venv/bin/python'
r = subprocess.run([PY, 'step_36.py'], cwd=REV, capture_output=True, text=True, timeout=3600)
print(r.stdout); print(r.stderr[-1000:] if r.stderr else "")
