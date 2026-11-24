
import subprocess, os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
REV = 'research/08_theory/D13r_lemma_review'
PY='<home>)/.venv/bin/python'
r = subprocess.run([PY, 'D13h_final_cell__step_10_J192.py'], cwd=WORK, capture_output=True, text=True, timeout=1200)
print("BUGGY FINAL CELL OUTPUT:")
print(r.stdout)
print("STDERR:", r.stderr[-500:] if r.stderr else "(none)")
