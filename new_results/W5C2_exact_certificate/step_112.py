import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 05:37:47 step_03 pm SOLVED: y_rat nz=11082/11101 (2023d/2012d), y_sq nz=11087 (2027d/2018d), 805s, cross-check OK.\n')
time.sleep(600)
print(subprocess.run(['bash','-c','grep -E "CONVERGED|y_rat|y_sq|cross-check|SAVED|step_03 mm|full recon" step_03.out | tail -12'], capture_output=True, text=True).stdout)