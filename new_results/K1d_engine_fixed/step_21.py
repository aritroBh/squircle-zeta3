s = open('step_42.py').read()
s = s.replace("assert 2**a2 * 3**a3 == rat, f\"log of {rat} not in {{log2,log3}}\"",
              "assert sp.Rational(2)**a2 * sp.Rational(3)**a3 == rat, f\"log of {rat} not in {{log2,log3}}\"")
open('step_42.py','w').write(s)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_42.py'],
                   capture_output=True, text=True, timeout=1200)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
