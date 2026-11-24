import subprocess
work = '<project>/research/08_theory/TOW3_tower_depth6'
for (delta, s) in [(3,0),(3,1),(4,0),(4,1),(4,2)]:
    r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','prove_adj56.py',str(delta),str(s)],
                       cwd=work, capture_output=True, text=True, timeout=1700)
    print(r.stdout.strip(), '| rc=', r.returncode)
