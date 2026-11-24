src = open("itab_build2.py").read()
old = """    out[m] = a.mid().str(dps=DPS+5, condense=DPS+10)"""
new = """    out[m] = a.mid().str(DPS+5)"""
assert old in src; src = src.replace(old, new)
old = "EPS = flint.arb(10)**(-(DPS+4)); STOP = (10.0)**(-(DPS-1))"
new = "EPS = flint.arb(10)**(-(DPS+15)); STOP = (10.0)**(-(DPS-1))"
assert old in src; src = src.replace(old, new)
open("itab_build2.py","w").write(src)
import subprocess, pickle, os
logf = open("itab_build.log", "w")
proc = subprocess.Popen(["<project>/.venv/bin/python", "-u", "itab_build2.py", "140"],
                        stdout=logf, stderr=logf, cwd=os.getcwd())
pickle.dump(proc.pid, open("itab_pid.pkl","wb"))
print("relaunched itab DPS=140, pid", proc.pid)
