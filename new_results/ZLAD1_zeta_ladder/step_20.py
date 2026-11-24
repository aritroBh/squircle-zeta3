import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_02.py','150'],
                   capture_output=True, text=True, timeout=1800)
print(r.stdout[-1200:])
print('ERR', r.stderr[-800:])
open('progress.log','a').write('step_02: exact mu_{k,n}, k<=6, n<=150 computed and validated vs G1(k<=3),G2(k=4) at n<=8; all rungs spot-checked\n')
