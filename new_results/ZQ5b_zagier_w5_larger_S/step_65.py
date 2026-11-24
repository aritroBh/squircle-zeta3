import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
s=open('step_09c_test.py').read()
s=s.replace("NS = flint.fmpz_mat(MT).nullspace()","NS, nulldim = flint.fmpz_mat(MT).nullspace()")
s=s.replace("n0 = NS.nrows(); n1 = NS.ncols()","n0 = NS.nrows(); n1 = NS.ncols(); print('nulldim', nulldim)")
open('step_09c_test.py','w').write(s)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','-u','step_09c_test.py'], capture_output=True, text=True, cwd='.', timeout=1200)
print(r.stdout); print("ERR:", r.stderr[-800:])
