import os, sys, subprocess, datetime
print("CWD:", os.getcwd())
print(os.listdir('.'))
os.makedirs('scratchpad', exist_ok=True)
# start progress log
with open('progress.log','a') as f:
    f.write(f"# AV8 run start {datetime.datetime.now().isoformat()}\n")
    f.write("step_00: initialized working directory\n")
print(sys.version)
import mpmath
print("mpmath", mpmath.__version__)
