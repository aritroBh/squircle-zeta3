code = '''"""step_06 (TOW5R): tau1/tau2(9),(10) from the DL definitions, exact Q(sqrt2,i) + 40-digit mpmath.
towlib defs: tau1(n) = (-1)^{n+1}[(i-b)^{-n-2}+(i+a)^{-n-2}-i^{-n-2}], a=1+r2, b=r2-1;
             tau2(n) = (-1)^{n+1}[(i-a)^{-n-2}+(i+b)^{-n-2}-i^{-n-2}].
Checks: exact QSI arithmetic -> Q(i) values == tau910.pkl; mpmath 40-digit; conj(tau1_n)=(-1)^n tau2_n;
        also n = 0..8 against the inherited table (sanity of the definition path).
"""
import os, sys, pickle
sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE + '/tow5_copy')
from fastred import QSI, IQ, AQ, BQ
from fractions import Fraction as Fr
from mpmath import mp, mpc, power
mp.mp.dps = 40

def tau1_qsi(n):
    s = n+1
    return ((-1)**s) * ((IQ-BQ)**(-(s+1)) + (IQ+AQ)**(-(s+1)) - IQ**(-(s+1)))
def tau2_qsi(n):
    s = n+1
    return ((-1)**s) * ((IQ-AQ)**(-(s+1)) + (IQ+BQ)**(-(s+1)) - IQ**(-(s+1)))

tau910 = pickle.load(open(_HERE + '/tow5_copy/tau910.pkl', 'rb'))
import sympy as sp
ok = True
msgs = []
for n in [9, 10]:
    t1 = tau1_qsi(n); t2 = tau2_qsi(n)
    # QSI = (p0 + p1 r2) + i (p2 + p3 r2); must land in Q(i): p1 = p3 = 0
    assert t1.p1 == 0 and t1.p3 == 0 and t2.p1 == 0 and t2.p3 == 0, f'sqrt2 survives at n={n}'
    v1 = sp.sympify(tau910[n][0], locals={'I': sp.I}); v2 = sp.sympify(tau910[n][1], locals={'I': sp.I})
    ok &= (Fr(sp.re(v1)) == t1.p0 and Fr(sp.im(v1)) == t1.p2)
    ok &= (Fr(sp.re(v2)) == t2.p0 and Fr(sp.im(v2)) == t2.p2)
    # conjugation symmetry
    ok &= (sp.expand(sp.conjugate(v1) - (-1)**n * v2) == 0)
    # 40-digit mpmath from the definition
    a = 1 + mp.sqrt(2); b = mp.sqrt(2) - 1; ii = mpc(0, 1)
    f1 = (-1)**(n+1) * (power(ii-b, -(n+2)) + power(ii+a, -(n+2)) - power(ii, -(n+2)))
    f2 = (-1)**(n+1) * (power(ii-a, -(n+2)) + power(ii+b, -(n+2)) - power(ii, -(n+2)))
    e1 = abs(f1 - (float(v1.as_real_imag()[0]) + 1j*float(v1.as_real_imag()[1])))
    e2 = abs(f2 - (float(v2.as_real_imag()[0]) + 1j*float(v2.as_real_imag()[1])))
    ok &= (e1 < mp.mpf('1e-35') and e2 < mp.mpf('1e-35'))
    msgs.append(f'tau1({n}) = {t1.p0} + ({t1.p2}) i; tau2({n}) = {t2.p0} + ({t2.p2}) i; mp residuals {mp.nstr(e1,3)},{mp.nstr(e2,3)}')
# sanity of the definition path at n = 0..8 vs inherited table tow2lib/tow4lib
import tow2lib as T2, tow4lib as T4
for n in range(0, 9):
    t1 = tau1_qsi(n); t2 = tau2_qsi(n)
    v1 = T4.tau1f(n); v2 = T4.tau2f(n)
    ok &= (t1.p1 == 0 and t1.p3 == 0 and Fr(sp.re(v1)) == t1.p0 and Fr(sp.im(v1)) == t1.p2)
    ok &= (t2.p1 == 0 and t2.p3 == 0 and Fr(sp.re(v2)) == t2.p0 and Fr(sp.im(v2)) == t2.p2)
print('STEP_06 VERDICT:', 'PASS' if ok else 'FAIL')
for m_ in msgs: print(' ', m_)
print('  conj(tau1_n) = (-1)^n tau2_n at n = 9, 10:', ok)
assert ok
'''
wd = '<project>/research/08_theory/TOW5R_review'
open(wd+'/step_06.py','w').write(code)
import subprocess, datetime
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'step_06.py'], cwd=wd, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print('STDERR tail:', r.stderr[-400:], 'rc:', r.returncode, 'elapsed:', datetime.datetime.now()-t0)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_06: tau(9),(10) exact QSI + 40-digit mp vs tau910.pkl + conj rules n<=10 rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
