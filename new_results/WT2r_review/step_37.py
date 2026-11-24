code = '''
import sys, time
sys.path.insert(0,'.')
from math import comb, factorial
from fractions import Fraction as Fr
from genintegrand import S_poly
import buildintegrand as bi
from fastred import *
IQ = QSI(0,0,1,0)
def sqrt_num_den(k, S):
    St, D = bi.Stilde(S)
    if not St: return None, None
    e = 2*k - 1 - D
    num = bi.qsi_pscale(St, bi.SQ2Q * QSI(Fr(-(8**k), 2*4**D)))
    dev = 0
    if e >= 0:
        for _ in range(e): num = bi.qsi_pmul(num, bi.VPOLY)
    else: dev = -e
    return num, dev
def qsi_peval(P, x):
    tot = ZERO
    for c in reversed(P): tot = tot*x + c
    return tot
def top_plus_i(k, S):
    num, dev = sqrt_num_den(k, S)
    if num is None: return ZERO
    ni = qsi_peval(num, IQ)
    d0 = (IQ**2 - AQ**2)*(IQ**2 - BQ**2) * (IQ**dev) * ((2*IQ)**(2*k-2))
    return ni / d0
for d in [10, 12, 14, 16]:
    t0=time.time()
    Ss = [S_poly(d,m) for m in range(d+1)]
    t1=time.time()
    tops = [top_plus_i(d, S) for S in Ss]
    t2=time.time()
    print(f"d={d}: S_polys {t1-t0:.1f}s, tops {t2-t1:.1f}s", flush=True)
'''
open('rev_prof5.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_prof5.py'], capture_output=True, text=True, cwd=work, timeout=300)
print(r.stdout); print(r.stderr[-500:] if r.returncode else "")
