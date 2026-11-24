src = '''"""step_74_local_check.py - Task 4: local block R fit on the 22-atom basis;
check z5 and z2z3 coefficients against prop:order5 at m=5.
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from q5_lll import fit
HERE = os.path.dirname(os.path.abspath(__file__))

def atoms(D):
    mp.mp.dps = D + 40
    phi = (1+mp.sqrt(5))/2
    l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
    L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
    L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
    K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
    K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
    K3m = L3m + lf**3/6
    A = {}
    for a in range(4):
        b = 3-a
        nm = 'z2' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        A[nm] = z2*l2**a*lf**b
    for a in range(3):
        b = 2-a
        nm = 'z3' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        A[nm] = z3*l2**a*lf**b
    A['z2^2l2'] = z2**2*l2; A['z2^2lf'] = z2**2*lf
    A['z5'] = z5; A['z2z3'] = z2*z3; A['z2L3m'] = z2*L3m
    A['K5m'] = K5m; A['K5q'] = K5q
    A['K4ml2'] = K4m*l2; A['K4mlf'] = K4m*lf; A['K4ql2'] = K4q*l2; A['K4qlf'] = K4q*lf
    A['K3ml2^2'] = K3m*l2**2; A['K3ml2lf'] = K3m*l2*lf; A['K3mlf^2'] = K3m*lf**2
    A['Lchi5'] = mp.mpf(pickle.load(open(os.path.join(HERE,'Lchi5_1550.pkl'),'rb'))['Lchi5_5'])
    INV = ['z2l2^3','z2l2^2lf','z2l2lf^2','z2lf^3','z3l2^2','z3l2lf','z3lf^2','z2^2l2','z2^2lf',
           'z5','z2z3','z2L3m','K5m','K5q','K4ml2','K4mlf','K4ql2','K4qlf','K3ml2^2','K3ml2lf','K3mlf^2','Lchi5']
    return A, INV

def main():
    mp.mp.dps = 800
    r5 = mp.sqrt(5); l2 = mp.log(2); lf = mp.log((1+r5)/2)
    W8 = pickle.load(open(os.path.join(HERE, 'wparts_dec_d800.pkl'),'rb'))
    lp = pickle.load(open(os.path.join(HERE, 'block_logparts.pkl'),'rb'))
    logv = sum((mp.mpf(p)+mp.mpf(q)*r5) * l2**a * lf**(5-a) for a, (p,q) in lp['local'].items())
    Rloc = mp.mpf(W8['blk5']['local']) - logv
    A, INV = atoms(800)
    mp.mp.dps = 800
    c = fit(Rloc, [(x, A[x]) for x in INV], cexp=int(800*0.92), nbest=2)
    h, res, rel, v0 = c[0]
    print("local R fit: height %.3e resid %s" % (h, mp.nstr(mp.mpf(res), 3)))
    print("z5 coefficient:", rel.get('z5'))
    print("z2z3 coefficient:", rel.get('z2z3'))
    print("atoms used:", sorted(rel.keys()))
    ok5 = rel.get('z5') == (__import__('fractions').Fraction(73125,64), __import__('fractions').Fraction(-30825,64))
    ok23 = rel.get('z2z3') == (__import__('fractions').Fraction(-64875,64), __import__('fractions').Fraction(32775,64))
    print("z5 == prop:order5 C_5^(5)/L = (73125/64, -30825/64):", ok5)
    print("z2z3 == prop:order5 C_5^(2,3)/L = (-64875/64, 32775/64):", ok23)
    print("TASK-4 CHECK:", "PASS" if (ok5 and ok23) else "FAIL")

if __name__ == "__main__":
    main()
'''
open('step_74_local_check.py','w').write(src)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'step_74_local_check.py'], capture_output=True, text=True, cwd='.')
print(r.stdout[-2000:])
if r.stderr: print("STDERR:", r.stderr[-400:])
print("wall %.0fs" % (time.time()-t0))
