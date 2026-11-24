work = '<project>/research/08_theory/TOW3R_review/'
import shutil
base = '<project>/research/08_theory/'
shutil.copy(base+'TOW2_tower_all_depths/genintegrand.py', work+'tow3_copy/genintegrand.py')
src = '''import sys, os, time, pickle
from fractions import Fraction as Fr
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work+'/tow3_copy')
import sympy as sp
from genintegrand import S_poly

v = sp.Symbol('v'); i = sp.I
aq = sp.sqrt(2)+1; bq = sp.sqrt(2)-1
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    den = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, den), gdiv(numB, den))

from collections import Counter
def parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from parts(rem-p, p, cur+[p])
def bell(delta, L):
    tot = sp.Integer(0)
    for part in parts(delta, 1, []):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz, mult in C.items():
            coef *= (L[sz]**mult) / sp.factorial(mult)
        tot += coef
    return tot

def tower_ratio(p, q, delta):
    return qsi_div(H[(p,q)][delta][p+q], H[(p,q)][0][p+q])

t0 = time.time()
for (p, q) in [(3,2),(5,2),(4,1),(3,3),(4,4),(7,3),(2,5),(6,6),(20,15),(0,5),(1,6)]:
    dd = p+q
    S = S_poly(dd, p)          # list of Fractions, low-to-high
    t = v/4 + sp.Rational(1,2) - 1/(4*v)
    St = sum(sp.Rational(c.numerator, c.denominator)*t**e for e, c in enumerate(S))
    logA = (2*dd-1)*sp.log(v) + sp.log(St) + (2-2*dd)*sp.log(v+i) - sp.log(v**2-aq**2) - sp.log(v**2-bq**2)
    L = {}
    cur = logA
    for r in range(1, 7):
        cur = sp.diff(cur, v)
        L[r] = sp.simplify(cur.subs({v: i})) / sp.factorial(r)
    ok = True
    detail = {}
    for delta in (5, 6):
        a_delta = sp.simplify(bell(delta, L))
        ratio = tower_ratio(p, q, delta)
        # ratio should be Gaussian rational = a_delta; sqrt2 part 0
        got = ratio[0]
        want_re = sp.re(a_delta); want_im = sp.im(a_delta)
        good = (ratio[1]==(Fr(0),Fr(0))) and (Fr(str(sp.nsimplify(want_re, rational=True)))==got[0] if not want_re.is_rational else sp.Rational(want_re)==got[0]) and (sp.Rational(want_im)==got[1])
        detail[delta] = good
        ok &= bool(good)
    print(f'key ({p},{q}) d={dd}: independent A(v) log-deriv -> a_5, a_6 match data ratios: {detail}  ({time.time()-t0:.0f}s)', flush=True)
'''
open(work+'step_06.py','w').write(src)
import subprocess, datetime
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,work+'step_06.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout, r.stderr[-2500:])
open(work+'step_06.out','w').write(r.stdout+r.stderr)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_06 independent A(v) log-deriv re-derivation at 11 keys\n")
