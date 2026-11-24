import subprocess, sys, os
print("out/d17:", sorted(os.listdir('out/d17')) if os.path.exists('out/d17') else "MISSING")
r = subprocess.run(['tail','-3','progress.log'],capture_output=True,text=True)
print(r.stdout)
