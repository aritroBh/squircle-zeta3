import subprocess, time
print(subprocess.run(["pgrep","-P","66142"],capture_output=True,text=True).stdout)
time.sleep(20)
r = subprocess.run(["pgrep","-P","66142"],capture_output=True,text=True).stdout
kids = [k for k in r.split() if k.strip()]
print("worker child PIDs of 66142:", kids, "count:", len(kids))
open("worker_pids.txt","w").write("\n".join(kids)+"\n")
