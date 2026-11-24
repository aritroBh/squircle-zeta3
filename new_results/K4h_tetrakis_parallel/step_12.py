import subprocess
r = subprocess.run(["pgrep","-fl","step_03.py"],capture_output=True,text=True)
print("before:", r.stdout)
subprocess.run(["pkill","-f","step_03.py"],capture_output=True,text=True)
import time; time.sleep(2)
r = subprocess.run(["pgrep","-fl","step_03.py"],capture_output=True,text=True)
print("after:", r.stdout or "(none)")
r = subprocess.run(["pgrep","-fl","python"],capture_output=True,text=True)
print("python procs:\n", r.stdout)
