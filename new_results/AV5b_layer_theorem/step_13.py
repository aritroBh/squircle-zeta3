import os, subprocess
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
print(os.listdir('.'))
print(open('step_05.py').read()[:200])
r = subprocess.run(["<project>/.venv/bin/python","step_05.py"],capture_output=True,text=True)
print(r.stdout)
print("ERR:", r.stderr[:2000])
