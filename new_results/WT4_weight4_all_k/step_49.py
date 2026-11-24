import subprocess
work='<project>/research/08_theory/WT4_weight4_all_k'
PY='<project>/.venv/bin/python'
for args in [['wt4_cert_step06.py','D2'],['wt4_cert_step06.py','D1']]:
    r=subprocess.run(['nice','-n','19',PY]+args, cwd=work, capture_output=True, text=True, timeout=1700)
    print('==', ' '.join(args))
    print(r.stdout)
    if r.returncode: print('STDERR:', r.stderr[-2500:]); break
