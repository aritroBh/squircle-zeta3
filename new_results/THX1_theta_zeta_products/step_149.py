
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', 'verify_thx1.py'],
                   capture_output=True, text=True, cwd='.')
lines = [l for l in r.stdout.split('\n') if l.startswith('PASS') or l.startswith('FAIL')]
print("PASS lines:", sum(1 for l in lines if l.startswith('PASS')))
print("FAIL lines:", sum(1 for l in lines if l.startswith('FAIL')))
import collections
groups = collections.Counter()
for l in lines:
    key = l.split(' ', 2)[1] if len(l.split(' ', 2)) > 1 else '?'
    groups[key.split(' ')[0]] += 1
print(groups)
print(r.stdout.strip().split('\n')[-1])
