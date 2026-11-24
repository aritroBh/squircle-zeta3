work = '<project>/research/08_theory/TOW3R_review/'
src = '''import sys, os, time
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work+'/tow3_copy'); os.chdir(work+'/tow3_copy')
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
w = sp.Symbol('w')
def qsqrt2_is_zero(expr):
    # expr rational in (d,m) over Q(i,sqrt2). Returns True iff expr == 0 in Q(i,sqrt2)(d,m).
    g = sp.together(expr)
    num = sp.fraction(g)[0]
    num = sp.expand(num.subs({sp.sqrt(2): w}))
    # reduce powers of w
    numr = sp.expand(num.subs({w: sp.sqrt(2)}))  # placeholder no-op
    # collect by parity of w-degree using series-free approach:
    even = sp.expand(num.subs({w: sp.Symbol('w2')}))  # not helpful; do it via Poly in w with extension disabled
    P = sp.Poly(num, w)
    a = sp.Integer(0); b = sp.Integer(0)
    for (e,), c in P.terms():
        if e % 2 == 0: a += c * 2**(e//2)
        else: b += c * 2**((e-1)//2)
    return sp.expand(a)==0 and sp.expand(b)==0, (a, b)
import pickle
L = pickle.load(open('L_parsed.pkl','rb'))
for r in range(1, 7):
    t0 = time.time()
    diff = L[r] - T.L(r)
    ok, (a,b) = qsqrt2_is_zero(diff)
    print(f'r={r}: L_parsed == tow2lib.L(r) in Q(i,sqrt2)(d,m): {ok}  ({time.time()-t0:.1f}s)', flush=True)
    if not ok:
        print('   even residual:', str(a)[:300]); print('   sqrt2 coeff:', str(b)[:300])
'''
open(work+'step_04.py','w').write(src)
import subprocess, datetime
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,work+'step_04.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout, r.stderr[-1500:])
open(work+'step_04.out','w').write(r.stdout+r.stderr)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_04 L equality Q(i,sqrt2) reduction r=1..6\n")
