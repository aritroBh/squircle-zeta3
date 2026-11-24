src=open('step_08_weight2_certificate.py').read()
old="""pp = []; start = 0
for i in range(r):"""
new="""pp = []; start = 0
for i in range(rk):"""
assert old in src
src=src.replace(old,new)
open('step_08_weight2_certificate.py','w').write(src)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_08_weight2_certificate.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout[-900:]); print('STDERR:', r.stderr[-800:])
open('progress.log','a').write('step_08: weight-2 CLOSED. New relation = Rogers five-term at x=sqrt2/2, y=sqrt2-1 (+4 more equivalent instances); quotient exactly 9 over Q; LLL hint spanned.\n')
