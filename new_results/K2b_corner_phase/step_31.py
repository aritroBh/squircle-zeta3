
wd = 'research/08_theory/K2b_corner_phase/'
code = open(wd+'local.py').read()
code = code.replace('ORD = 5','ORD = 8')
extra = r'''
g3 = sp.expand(Bser.coeff(sq,7)/gam)
g4 = sp.expand(Bser.coeff(sq,9)/gam)
print("g3 =", sp.simplify(g3.subs(Eb,1/E)))
print("g4 =", sp.simplify(g4.subs(Eb,1/E)))
'''
open(wd+'local2.py','w').write(code + extra)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python', wd+'local2.py'], capture_output=True, text=True, timeout=900)
print(r.stdout[-1500:]); print(r.stderr[-1500:] if r.stderr else '')
