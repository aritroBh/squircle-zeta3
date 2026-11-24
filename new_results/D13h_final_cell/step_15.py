src = open('step_10.py').read()
old = """assert (al_ + 2*L_J) > 6, "monotonicity windows for y*w^q fail\""""
new = """assert (al_ + 2*L_J) > 6, "monotonicity windows for y*w^q fail"
sqyJ0 = y_J.sqrt()   # sqrt(y_J) point value (box sqrt of [0,y_J] is nan in arb)"""
assert old in src; src = src.replace(old, new)
old = """sqyJ0 = y_J.sqrt()
sqyJ = sqyJ0"""
new = """sqyJ = sqyJ0"""
assert old in src; src = src.replace(old, new)
open('step_10.py','w').write(src)
import subprocess, sys
r = subprocess.run([sys.executable, 'step_10.py'], capture_output=True, text=True, timeout=300)
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
