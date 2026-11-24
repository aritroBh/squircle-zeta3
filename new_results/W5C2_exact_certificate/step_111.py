import subprocess
print(subprocess.run(['bash','-c','grep -E "CONVERGED|y_rat|y_sq|cross-check|SAVED|step_03 (pm|mp|mm):" step_03.out | tail -20'], capture_output=True, text=True).stdout)