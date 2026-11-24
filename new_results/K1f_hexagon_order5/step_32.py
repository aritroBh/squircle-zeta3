import subprocess
PY='<project>/.venv/bin/python'
wd = '<project>/research/08_theory/K1f_hexagon_order5'
r = subprocess.run([PY, f'{wd}/step_18.py'], capture_output=True, text=True, cwd=wd)
print(r.stdout)
print(r.stderr[-2000:] if r.returncode else '')
