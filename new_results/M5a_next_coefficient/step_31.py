import subprocess, time
work = '<project>/research/08_theory/M5a_next_coefficient'
PY = '<project>/.venv/bin/python'
t0 = time.time()
for target in ['w6', 'm5']:
    for dps in ['520', '470']:
        r = subprocess.run(['nice', PY, '-B', 'step_hp5.py', target, dps], capture_output=True, text=True, cwd=work)
        print(r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr[:500])
print("total elapsed:", time.time()-t0)
with open(f'{work}/progress.log', 'a') as f:
    f.write("step: step_hp5.py w6/m5 at working 520 and 470 dps -> hp5_results.json\n")