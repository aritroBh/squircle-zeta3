import os, subprocess
cwd = os.getcwd()
# check background values progress
pid = open('step_12_pid.txt').read().strip()
r = subprocess.run(['ps','-p',pid,'-o','pid,etime,%cpu,comm'], capture_output=True, text=True)
print("PS:", r.stdout)
if os.path.exists('step_12_bg.log'):
    print("BG LOG (tail):")
    print(open('step_12_bg.log').read()[-1500:])
