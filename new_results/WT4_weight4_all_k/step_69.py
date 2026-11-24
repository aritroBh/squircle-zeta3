import subprocess, os, json
work='<project>/research/08_theory/WT4_weight4_all_k'
PY='<project>/.venv/bin/python'
# capture verify_wt4.py output for the record
r=subprocess.run(['nice','-n','19',PY,'verify_wt4.py'], cwd=work, capture_output=True, text=True, timeout=1700)
open(os.path.join(work,'verify_wt4.out'),'w').write(r.stdout)
print(r.stdout[-600:])
print('returncode', r.returncode)
