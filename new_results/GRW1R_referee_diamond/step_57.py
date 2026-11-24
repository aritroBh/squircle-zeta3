import subprocess, time
t0=time.time()
p = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','verify_grw1r.py'],
                   capture_output=True, text=True, timeout=1750)
tail = p.stdout.split('==== [6]')[-1]
print('...'+tail[:1500])
print("STDERR tail:", p.stderr[-500:])
print("wall:", time.time()-t0)
with open('progress.log','a') as f:
    f.write("step_35: final verify_grw1r.py end-to-end run OK\n")
