
import subprocess, os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
for f in ['verify_Btrue.py', 'certify_H2.py']:
    r = subprocess.run([PY, f], cwd=FIX, capture_output=True, text=True, timeout=1800)
    print("="*20, f, "="*20)
    print(r.stdout); print(r.stderr[-800:] if r.stderr else "")
