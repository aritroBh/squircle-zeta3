import subprocess, time
t0=time.time()
p = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','verify_grw1r.py'],
                   capture_output=True, text=True, timeout=1750)
print(p.stdout)
print("STDERR:", p.stderr[-2000:])
print("wall:", time.time()-t0)
with open('progress.log','a') as f:
    f.write(f"step_31: verify_grw1r.py run OK\n")
