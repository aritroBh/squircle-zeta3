step18 = r'''"""step_18.py -- fit T(r) in a candidate span via overdetermined least squares + holdout check."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 60

def T_alpha_c(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def basis(c):
    r = (1-c)/(1+c)
    L = mp.log; A = lambda x: mp.polylog(2, x)
    return [('1', mp.mpf(1)), ('pi^2', mp.pi**2),
            ('Li2(r)', A(r)), ('Li2(-r)', A(-r)), ('Li2(c)', A(c)), ('Li2(-c)', A(-c)),
            ('Li2((1+r)/2', A((1+r)/2)), ('Li2((1-r)/2', A((1-r)/2)),
            ('Li2(2r/(1+r)', A(2*r/(1+r))), ('Li2((1-c)/2', A((1-c)/2)), ('Li2((1+c)/2', A((1+c)/2)),
            ('Li2(r/(1+r)', A(r/(1+r))), ('Li2(-r/(1+r)', A(-r/(1+r))),
            ('log2^2', L(2)**2), ('logr^2', L(r)**2), ('log(1+r)^2', L(1+r)**2), ('log(1-r)^2', L(1-r)**2),
            ('logc^2', L(c)**2), ('log(1+c)^2', L(1+c)**2), ('log(1-c)^2', L(1-c)**2),
            ('logrlog(1+r)', L(r)*L(1+r)), ('logrlog(1-r)', L(r)*L(1-r)), ('logrlog2', L(r)*L(2)),
            ('logclog(1+c)', L(c)*L(1+c)), ('logclog2', L(c)*L(2)), ('log(1+r)log2', L(1+r)*L(2))]

fitangles = [mp.mpf(x) for x in (mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7,
               mp.mpf(5)/6, mp.mpf(1)/5, mp.mpf(3)/7, mp.mpf(5)/8, mp.mpf(2)/9, mp.mpf(5)/7, mp.mpf(7)/9,
               mp.mpf(1)/6, mp.mpf(3)/8, mp.mpf(4)/7, mp.mpf(5)/9, mp.mpf(6)/7, mp.mpf(7)/8)]
holdout = [mp.mpf(1)/7, mp.mpf(4)/9, mp.mpf(2)/3, mp.mpf(6)/11]
print("computing T at", len(fitangles), "angles...")
Tfit = [T_alpha_c(c) for c in fitangles]
Thold = [T_alpha_c(c) for c in holdout]
Bfit = [[x for _, x in basis(c)] for c in fitangles]
Bhold = [[x for _, x in basis(c)] for c in holdout]
names = [n for n, _ in basis(fitangles[0])]
N = len(names)
A = mp.matrix(Bfit); Tv = mp.matrix(Tfit)
# least squares via normal equations with SVD-like: use mp.lu_solve on A^T A (may be ill-conditioned); use qr
import mpmath
# use Householder: mp doesn't expose; do SVD via numpy-free: use mp.svd? use matrix ops
AtA = A.T * A
AtT = A.T * Tv
try:
    q = mp.lu_solve(AtA, AtT)
except Exception as ex:
    print("lu fail", ex); q = None
if q is not None:
    resid = (A*q - Tv)
    print("LS residual max:", mp.nstr(mp.matrix(list(resid)).norm('inf'), 3))
    # rationalize
    qq = []
    for x in q:
        f = mp.fraction(x, 1) if False else None
        qq.append(x)
    print("coefficients (float):")
    for n, x in zip(names, q):
        if abs(x) > mp.mpf('1e-40'):
            print(f"   {n}: {mp.nstr(x, 30)}")
print("STEP18 DONE")
'''
open('step_18.py','w').write(step18)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_18.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
