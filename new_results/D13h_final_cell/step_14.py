src = open('step_10.py').read()
old = """fB = cf(ERA, Y, R); gB = cg(ERA, Y, R); fpB = cfp(ERA, Y, R); gpB = cgp(ERA, Y, R)
Lip_eta = up(pi*rU**2/4)*(up(fpB)*up(gB) + up(fB)*up(gpB))/lo(gB)**2"""
new = """fB = cf(ERA, Y, R); gB = cg(ERA, Y, R); fpB = cfp(ERA, Y, R); gpB = cgp(ERA, Y, R)
diffB = fpB*gB - fB*gpB            # = g^2 d/deta(f/g), NEGATIVE (curv2reg decreasing in eta)
Lip_eta = up(pi*rU**2/4)*max(abs(lo(diffB)), abs(up(diffB)))/lo(gB)**2"""
assert old in src; src = src.replace(old, new)
old = "E1R = (Lip_eta*up(Y.sqrt())*w_J**3*rU/(2*c2reg_lo)**arb(\"1.5\")).exp()"
new = "E1R = (Lip_eta*sqyJ0*w_J**3*rU/(2*c2reg_lo)**arb(\"1.5\")).exp()"
assert old in src; src = src.replace(old, new)
old = "sqyJ = up(Y.sqrt())"
new = "sqyJ0 = y_J.sqrt()\nsqyJ = sqyJ0"
assert old in src; src = src.replace(old, new)
open('step_10.py','w').write(src)
import subprocess, sys
r = subprocess.run([sys.executable, 'step_10.py'], capture_output=True, text=True, timeout=300)
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
