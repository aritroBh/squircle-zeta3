import subprocess, os
work='<project>/research/08_theory/WT4_weight4_all_k'
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','19',PY,'wt4_A_kernels.py'], cwd=work, capture_output=True, text=True, timeout=1500)
print(r.stdout)
print('STDERR:', r.stderr[-2000:])
