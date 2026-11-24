import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 05:50:45 step_03 mp SOLVED: y_rat nz=11037/11041 (1991d/1981d), y_sq nz=11038 (1991d/1982d), 776s, cross-check OK.\n')
time.sleep(600)
print(subprocess.run(['bash','-c','grep -E "step_03 mm|CONVERGED|y_rat|y_sq|cross-check|SAVED|DONE" step_03.out | tail -12'], capture_output=True, text=True).stdout)