step34 = r'''import sys, json, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
from cor13c_common import cross_table
from fractions import Fraction as Fr
tab = cross_table()
vals = sorted(tab.keys(), key=lambda z: (z.denominator, z.numerator))
print("cross-ratio values of {0,1/3,1,-1,inf}:", [str(z) for z in vals])
# group into S3 orbits
def orbit(l):
    return {l, 1-l, 1/l, 1/(1-l), l/(l-1), (l-1)/l}
seen = set(); orbits = []
for z in vals:
    if z in seen: continue
    o = orbit(z); seen |= o; orbits.append(sorted(o, key=lambda w:(w.denominator,w.numerator)))
print("S3 orbits:", [[str(z) for z in o] for o in orbits])
# certificate stats
import os
for f in ['cert13_w23.json', 'cert13_w4.json', 'cert13_euler.json']:
    print(f, os.path.getsize(f), 'bytes')
c23 = json.load(open('cert13_w23.json'))
for n in ('2','3'):
    ys = [Fr(int(a), int(b)) for a, b in c23[n]['y']]
    hn = [max(len(str(abs(c.numerator))), len(str(c.denominator))) for c in ys]
    nz = sum(1 for c in ys if c)
    print(f"w{n}: rows={len(c23[n]['prow'])} multipliers={len(ys)} nonzero={nz} max_height={max(hn)} digits")
y4 = pickle.load(open('w4_y_13.pkl','rb'))
hn = [max(len(str(abs(c.numerator))), len(str(c.denominator))) for c in y4]
nz = sum(1 for c in y4 if c)
print(f"w4: rows=1695 multipliers={len(y4)} nonzero={nz} max_height={max(hn)} digits median={sorted(hn)[len(hn)//2]}")
with open('progress.log','a') as f:
    f.write('step_34: report data gathered (orbits, cert sizes, heights)\n')
'''
open('step_34.py','w').write(step34)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', 'step_34.py'], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-2000:] if r.stderr else "")
