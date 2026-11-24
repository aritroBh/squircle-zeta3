import subprocess, time
time.sleep(420)
print(subprocess.run(['bash','-c','grep -E "pmb|step_01b pm|mpb|step_01b mp|plateau|DONE|skip" step_01b.out | tail -12'], capture_output=True, text=True).stdout)