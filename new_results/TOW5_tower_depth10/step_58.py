import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
for (dl, s) in [(9,2),(9,3),(9,4)]:
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY,'prove_adj910_fast.py',str(dl),str(s)], cwd=work, capture_output=True, text=True, timeout=1700)
    print(r.stdout.strip())
    if r.stderr: print("ERR", r.stderr[-300:])
    print(f"  -> ({dl},{s}) elapsed {round(time.time()-t0,1)}s rc={r.returncode}")
    with open(work+'/progress.log','a') as f:
        f.write(f"{time.strftime('%H:%M:%S')} step: prove_adj910_fast.py {dl} {s} rc={r.returncode} ({round(time.time()-t0,1)}s)\n")
