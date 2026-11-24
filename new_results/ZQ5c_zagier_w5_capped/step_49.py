import os, time, subprocess
cwd=os.getcwd()
# poll until done (check for 'DONE' in log or pkl exists with cross-check)
for _ in range(20):
    log = open('step_19_bg.log').read()
    if 'D13 VALS DONE' in log:
        break
    time.sleep(30)
print(log)
