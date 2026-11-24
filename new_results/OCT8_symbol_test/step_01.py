import os, subprocess, textwrap, json, time
work = "<project>/OCT8_symbol"
os.makedirs(work, exist_ok=True)
os.chdir(work)
# check siblings
print(os.listdir(".."))
with open("progress.log","a") as f:
    f.write(f"=== OCT8 session start {time.ctime()} ===\n")
PY='<project>/.venv/bin/python'
r = subprocess.run([PY,"-c","import mpmath, sympy; print(mpmath.__version__, sympy.__version__)"],capture_output=True,text=True)
print(r.stdout, r.stderr)
