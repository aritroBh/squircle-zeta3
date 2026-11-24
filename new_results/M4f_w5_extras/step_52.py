import subprocess, time
WORK='<project>/research/08_theory/M4f_w5_extras'
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','verify_json.py','R5_pareto.json','R5_min3.json'],
                   cwd=WORK, capture_output=True, text=True, timeout=1200)
lines = [l for l in r.stdout.splitlines() if l.startswith(('1.','2a.','2b.','3.'))]
print('verify pareto+min3:', *lines, sep='\n')
print('STDERR', r.stderr[-300:] if r.stderr else '(clean)')
t1=time.time()
r2 = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','overcomplete.py','1500'],
                   cwd=WORK, capture_output=True, text=True, timeout=2400)
print('1500 run elapsed', time.time()-t1)
print(r2.stdout[-2500:])
print('STDERR', r2.stderr[-500:])
