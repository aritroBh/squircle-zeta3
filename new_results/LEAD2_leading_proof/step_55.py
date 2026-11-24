import subprocess, os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'verify_lead2.py'], capture_output=True, text=True, timeout=900)
print(r.stdout[-3000:])
print(f"elapsed {time.time()-t0:.0f}s rc={r.returncode}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: final full verify_lead2.py run (A-G, d<=11) -> rc={r.returncode}, {time.time()-t0:.0f}s\n")
