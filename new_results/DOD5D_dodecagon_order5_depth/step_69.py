import subprocess, pickle, os
pid = pickle.load(open("itab_pid.pkl","rb"))
try: os.kill(pid, 9)
except Exception as e: print("kill:", e)
src = open("itab_build2.py").read()
old = "EPS = flint.arb(10)**(-(DPS-8)); STOP = (10.0)**(-(DPS-5))"
new = "EPS = flint.arb(10)**(-(DPS+4)); STOP = (10.0)**(-(DPS-1))"
assert old in src; src = src.replace(old, new)
open("itab_build2.py","w").write(src)
logf = open("itab_build.log", "w")
proc = subprocess.Popen(["<project>/.venv/bin/python", "-u", "itab_build2.py", "135"],
                        stdout=logf, stderr=logf, cwd=os.getcwd())
pickle.dump(proc.pid, open("itab_pid.pkl","wb"))
print("relaunched itab at DPS=135, pid", proc.pid)
print("--- gtab log ---")
print(open("gtab_build.log").read())
