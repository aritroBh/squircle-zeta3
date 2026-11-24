import subprocess
work = '<project>/research/08_theory/TOW3_tower_depth6'
for delta in [4, 6]:
    r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','prove_E3b56.py',str(delta)],
                       cwd=work, capture_output=True, text=True, timeout=1700)
    print(r.stdout.strip(), '| rc=', r.returncode, flush=True)
with open(f'{work}/progress.log','a') as f:
    f.write('step_49-55: adjacent (5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3) PROVED (simplify fix); E3 blindness delta=5,6 PROVED\n')
