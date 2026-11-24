s = open('step_43.py').read()
s = s.replace("""    from fractions import Fraction
    prop = all(Fraction(a,b)==Fraction(reln[0],exactn[0]) for a,b in zip(reln,exactn)) if exactn[0]!=0 else False""",
"""    prop = all(a*exactn[0]==b*reln[0] for a,b in zip(reln,exactn))""")
open('step_43.py','w').write(s)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_43.py'],
                   capture_output=True, text=True, timeout=900)
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "")
