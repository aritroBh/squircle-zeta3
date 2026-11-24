import subprocess
WORK='<project>/research/08_theory/M4f_w5_extras'
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','drop_L.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=600)
print(r.stdout)
print('STDERR', r.stderr[-500:])
with open(f'{WORK}/progress.log','a') as fh:
    fh.write('verify_json: all 3 forms PASS (symbol exact, cobracket exact, -2.1052657e-1596 each)\n')
    fh.write('overcomplete 1500: same 14-dep span (combined rank 14), identical target relation, resid 2.17e-1589\n')
