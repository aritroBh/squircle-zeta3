src=open('step_08_weight2_certificate.py').read()
old = """M = mat(known + [new_row]); rr, r = M.rref(inplace=True)
# membership test: add lll row, rank must stay r"""
new = """M = mat(known + [new_row]); _, r = M.rref(inplace=True)
Mk = mat(known); rr, rk = Mk.rref(inplace=True)
assert rk == r - 1
# membership test: add lll row, rank must stay r"""
assert old in src
src = src.replace(old, new)
open('step_08_weight2_certificate.py','w').write(src)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_08_weight2_certificate.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout[-1400:]); print('STDERR:', r.stderr[-1200:])
