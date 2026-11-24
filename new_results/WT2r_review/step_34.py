code = '''
import sys, time
sys.path.insert(0,'.')
from genintegrand import S_poly
for d in [10,12,14,16,18,20]:
    t0=time.time()
    Ss = [S_poly(d,m) for m in range(0,d+1,5)]
    print(f"d={d}: 5 S_polys in {time.time()-t0:.1f}s")
'''
open('rev_prof2.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_prof2.py'], capture_output=True, text=True, cwd=work, timeout=600)
print(r.stdout); print(r.stderr[-500:] if r.returncode else "")
