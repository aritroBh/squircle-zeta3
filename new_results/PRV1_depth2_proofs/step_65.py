import subprocess, time
PY = '<project>/.venv/bin/python'
wd = 'research/08_theory/PRV1_depth2_proofs'

# 1. CX1's original verifier (the input)
cx = 'research/08_theory/CX1_octagon_depth2_proof/verify.py'
t0 = time.time()
r1 = subprocess.run([PY, cx], capture_output=True, text=True, timeout=1800)
print("CX1 verify.py tail:", r1.stdout.strip().split('\n')[-3:], f"({time.time()-t0:.0f}s)")
print("CX1 stderr:", r1.stderr[-300:])

# 2. Final clean run of verify_full.py
t0 = time.time()
r2 = subprocess.run([PY, f'{wd}/verify_full.py'], capture_output=True, text=True, cwd=wd, timeout=1800)
print("verify_full.py tail:", r2.stdout.strip().split('\n')[-2:], f"({time.time()-t0:.0f}s)")
print("stderr:", r2.stderr[-300:])
with open(f'{wd}/verify_full.out','w') as f:
    f.write(r2.stdout)
ok1 = 'ALL CHECKS PASSED' in r1.stdout
ok2 = 'ALL CHECKS PASSED' in r2.stdout
print("CX1 PASSED:", ok1, "| verify_full PASSED:", ok2)
