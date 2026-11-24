import subprocess
r = subprocess.run(["python3","q1s_exact.py"],capture_output=True,text=True)
print("OUT:", r.stdout); print("ERR:", r.stderr[:2000])
