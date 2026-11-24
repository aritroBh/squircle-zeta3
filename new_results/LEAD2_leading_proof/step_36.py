import subprocess, os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'verify_lead2.py'], capture_output=True, text=True, timeout=900)
print(r.stdout)
print('STDERR tail:', r.stderr[-800:])
print(f"elapsed {time.time()-t0:.0f}s rc={r.returncode}")
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: reran verify_lead2.py after tol fixes -> rc={r.returncode}\n")
