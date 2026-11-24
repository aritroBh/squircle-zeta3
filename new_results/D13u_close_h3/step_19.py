import subprocess, os, json, math
r = subprocess.run(['<project>/.venv/bin/python',
    '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed/step_18.py'],
    capture_output=True, text=True, cwd='/tmp')
print("verbatim step_18.py stdout:"); print(r.stdout)
print("stderr:", r.stderr[-300:] if r.stderr else "(none)")
# compare with stored json
stored = json.load(open('<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed/tight_100_400.json'))
mine = json.load(open('/tmp/tight_100_400.json'))
dmax = max(abs(m[2]-s[2]) for m,s in zip(mine, stored))
print(f"max |C diff| vs stored tight_100_400.json over {len(stored)} cells: {dmax:.3e}")
