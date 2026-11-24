code = '''
import sys, time
sys.path.insert(0,'.')
from math import comb, factorial
from fractions import Fraction as Fr
from genintegrand import S_poly
import buildintegrand as bi
from fastred import *
IQ = QSI(0,0,1,0)
t0=time.time()
for m in range(0, 11):
    tm=time.time()
    S = S_poly(10,m)
    St, D = bi.Stilde(S)
    print(f"m={m}: S_poly+Stilde {time.time()-tm:.2f}s (D={D}, len(St)={len(St)})", flush=True)
print("total", time.time()-t0)
'''
open('rev_prof4.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_prof4.py'], capture_output=True, text=True, cwd=work, timeout=120)
print(r.stdout); print(r.stderr[-500:] if r.returncode else "")
