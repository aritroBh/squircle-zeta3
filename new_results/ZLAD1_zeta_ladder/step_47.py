import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r2 = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','verify_zlad1.py'],
                   capture_output=True, text=True, timeout=600)
print(r2.stdout)
print('exit', r2.returncode)
