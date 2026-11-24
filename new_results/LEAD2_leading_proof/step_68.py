import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
i4 = txt.find('\n## 4.'); i7 = txt.find('\n## 7.')
open('proof_part1.md','w').write(txt[:i4])
open('proof_part2.md','w').write(txt[i4:i7])
assert txt[:i4] + txt[i4:i7] + txt[i7:] == txt
print('part files synced exactly; master length', len(txt))
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'verify_lead2.py'], capture_output=True, text=True, timeout=900)
print(r.stdout.strip().split('\n')[-1])
print('rc=',r.returncode)
with open('progress.log','a') as f:
    import time
    f.write(f"{time.strftime('%H:%M:%S')} step: part files synced; verify_lead2.py rerun rc={r.returncode}\n")
