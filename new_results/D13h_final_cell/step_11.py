# step_11.py - fix the d<x0/3 assert (avoid log of box containing 0; use endpoint monotonicity)
src = open('step_10.py').read()
old = """# d = w/sqrt(2 c2) = w / sqrt(2 y c2reg)  [assert d < x0/3 i.e. y w^2 < 2 c2reg/(9 rho^2)]
assert up(Y*(al_ + 2*(-Y.log()))) < 2*c2reg_lo/(9*rU**2), "d < x0/3 fails\""""
new = """# d = w/sqrt(2 c2) = w / sqrt(2 y c2reg)  [assert d < x0/3 i.e. y w^2 < 2 c2reg/(9 rho^2)]
# y w^2 = y(al_+2L) is increasing on (0,y_J] since al_+2L_J > 2 -> bound by endpoint:
yw2_up = y_J*(al_ + 2*L_J)
assert yw2_up < 2*c2reg_lo/(9*rU**2), "d < x0/3 fails\""""
assert old in src
src = src.replace(old, new)
open('step_10.py','w').write(src)
import subprocess, sys
r = subprocess.run([sys.executable, 'step_10.py'], capture_output=True, text=True, timeout=300)
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
