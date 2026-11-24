import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_01.py'], capture_output=True, text=True, cwd=work)
print(r.stdout)
print(r.stderr[-500:] if r.returncode else "")
