import subprocess
print(subprocess.run(['bash','-c','sed -n 1,120p "<project>/research/08_theory/E5c_nonlocal_z5/mzv.py"'], capture_output=True, text=True).stdout)
