import subprocess
work='<project>/research/08_theory/WT4_weight4_all_k'
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','19',PY,'verify_wt4.py'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR:', r.stderr[-1500:])
print('returncode', r.returncode)
