import os, subprocess
print("PATH:", os.environ.get("PATH"))
print("which python3:", subprocess.run(["which","python3"],capture_output=True,text=True).stdout)
print("VIRTUAL_ENV:", os.environ.get("VIRTUAL_ENV"))
import sys; print("tool python:", sys.executable)
r=subprocess.run([sys.executable,"-c","import mpmath;print(mpmath.__version__)"],capture_output=True,text=True)
print("venv mpmath:", r.stdout, r.stderr)
