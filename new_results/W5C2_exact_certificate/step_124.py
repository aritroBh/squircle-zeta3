import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step06 = r'''"""step_06: is every Galois transform g.r (g in {tau,sigma,sigma*tau}) of every used unfolded row
literally present (as a vector, up to nonzero scalar) in the unfolded row files? If yes for all,
then D = explicit combination of file rows (no covariance argument needed)."""
import sys, os, time, pickle, json
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from math import gcd
from w5c2_fold import Ctx, gapply, W5C1

def log(s):
    s = time.strftime('%Y-%m-%d %H:%M:%S ')+s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

SECNAME = {(1,1):'pp',(1,-1):'pm',(-1,1):'mp',(-1,-1):'mm'}
t0 = time.monotonic()
ctx = Ctx(); IDS = ctx.IDS
def normt(pairs):
    d = {}
    for w, a in pairs:
        d[w] = d.get(w, 0) + a
    d = {w: a for w, a in d.items() if a}
    if not d: return ()
    g = 0
    for a in d.values(): g = gcd(g, abs(a))
    items = sorted(d.items())
    if items[0][1] < 0: g = -g
    return tuple((w, a//g) for w, a in items)
ENG = ['rels_std.pkl','rels_lift1.pkl','rels_lift2old.pkl','rels_lift3old.pkl','rels_lift4old.pkl','rels_frame234.pkl']
rowset = set()
for fn in ENG:
    R = pickle.load(open(os.path.join(W5C1, fn),'rb'))
    for lab, row in R:
        rowset.add(normt(row))
VW5 = pickle.load(open(os.path.join(W5C1,'validated_w5.pkl'),'rb'))['rows']
for row in VW5:
    rowset.add(normt([(IDS[w], c) for w, c in row]))
log(f'step_06: unfolded row set built {len(rowset)} unique normalized rows ({time.monotonic()-t0:.0f}s)')
# used unique unfolded rows
used = {}
for sec, tag in SECNAME.items():
    cert = json.load(open(f'certificate_{tag}.json'))
    for ent in cert['rows']:
        ufam, uidx = ent['src']['ufam'], ent['src']['uidx']
        used[(ufam, uidx)] = True
log(f'step_06: unique used unfolded rows {len(used)}')
missing = []
present = 0
UROWS = {fam: pickle.load(open(os.path.join(W5C1, fn),'rb')) for fam, fn in
         {'std':'rels_std.pkl','lift1':'rels_lift1.pkl','lift2old':'rels_lift2old.pkl',
          'lift3old':'rels_lift3old.pkl','lift4old':'rels_lift4old.pkl','frame234':'rels_frame234.pkl'}.items()}
for (ufam, uidx) in used:
    if ufam == 'framew5':
        urow = [(IDS[w], c) for w, c in VW5[uidx]]
    else:
        urow = UROWS[ufam][uidx][1]
    for g in (1, 2, 3):
        t = normt([(gapply(w, g), a) for w, a in urow] if False else [(IDS[gapply(tuple(ctx.W5[w]), g)], a) for w, a in urow])
        if t in rowset: present += 1
        else: missing.append((ufam, uidx, g))
log(f'step_06: g.r present in row files: {present}/{3*len(used)}; missing={len(missing)} ({time.monotonic()-t0:.0f}s)')
from collections import Counter
mc = Counter((m[0], m[2]) for m in missing)
log(f'step_06: missing by (family, g): {dict(mc)}')
'''
open('step_06.py','w').write(step06)
import subprocess, time
out = open('step_06.out','wb')
p = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','step_06.py'],
                     cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
open('step_06.pid','w').write(str(p.pid))
print('launched step_06, pid', p.pid)
time.sleep(120)
print(open('step_06.out').read())