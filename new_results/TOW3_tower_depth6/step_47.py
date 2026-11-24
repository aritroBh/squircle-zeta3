import subprocess, time
work = '<project>/research/08_theory/TOW3_tower_depth6'
for (delta, s) in [(6,1),(6,2),(6,3)]:
    r = subprocess.run(['nice','-n','10', '<project>/.venv/bin/python', 'prove_V56.py', str(delta), str(s)],
                       cwd=work, capture_output=True, text=True, timeout=1700)
    print(r.stdout.strip(), r.stderr[-300:] if r.returncode else '', flush=True)
with open(f'{work}/progress.log','a') as f:
    f.write('step_44-47: V(5),V(6) + slopes (5,1),(6,1),(6,2) + tops (5,2),(6,3) non-adjacent PROVED (numerator 0)\n')
