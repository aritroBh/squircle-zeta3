import subprocess, time, os
t0 = time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_cor13c.py'], capture_output=True, text=True)
el = time.time()-t0
lines = r.stdout.strip().split('\n')
print('\n'.join(lines[:8]))
print('...')
print('\n'.join(lines[-6:]))
fails = [l for l in lines if 'FAIL' in l]
passes = sum(1 for l in lines if '[PASS]' in l)
print(f"\nfresh run: {el:.1f}s, {passes} PASS, {len(fails)} FAIL, exit={r.returncode}")
with open('progress.log','a') as f:
    f.write(f'step_35: fresh verify_cor13c.py run {el:.0f}s, {passes} PASS, {len(fails)} FAIL\n')
# also list files for the report
print("\nfiles:")
for p in sorted(os.listdir('.')):
    if p.endswith(('.py','.json','.log','.pkl','.txt')) or p=='result.md':
        print(f"  {p} {os.path.getsize(p)}")
