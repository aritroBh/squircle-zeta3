code = '''
import sys, time
sys.path.insert(0,'.')
from math import comb, factorial
import buildintegrand
from genintegrand import S_poly
from buildintegrand import sqrt_part_PF
from fastred import *

t0=time.time()
for d in range(2,13):
    td=time.time()
    Ps = [sqrt_part_PF(d,S_poly(d,m)) for m in range(d+1)]
    print(f"d={d}: built {len(Ps)} sqrt parts in {time.time()-td:.1f}s (total {time.time()-t0:.1f}s), max ipos len {max(len(P[1]) for P in Ps)}")
'''
open('rev_prof.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_prof.py'], capture_output=True, text=True, cwd=work, timeout=300)
print(r.stdout); print(r.stderr[-500:] if r.returncode else "")
