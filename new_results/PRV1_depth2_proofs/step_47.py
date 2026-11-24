import subprocess
PY = '<project>/.venv/bin/python'
wd = 'research/08_theory/PRV1_depth2_proofs'
r = subprocess.run([PY, f'{wd}/_trial.py'], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout[-5500:])
print("STDERR:", r.stderr[-1500:])
