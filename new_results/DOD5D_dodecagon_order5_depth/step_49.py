import subprocess, os, time
# kill any stragglers, then launch in background
subprocess.run(["pkill", "-f", "itab_build"], capture_output=True)
logf = open("itab_build.log", "w")
proc = subprocess.Popen(["<project>/.venv/bin/python", "-u", "itab_build2.py", "130"],
                        stdout=logf, stderr=logf, cwd=os.getcwd())
pickle = __import__("pickle")
pickle.dump(proc.pid, open("itab_pid.pkl","wb"))
print("launched pid", proc.pid)
