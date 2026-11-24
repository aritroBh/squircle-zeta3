step20 = r'''"""step_20.py -- simultaneous identification of T(r): greedy independent subset + exact solve + holdout."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 100

def T_alpha_c(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def cand(c):
    r = (1-c)/(1+c)
    L = mp.log; A = lambda x: mp.polylog(2, x)
    return [('pi^2', mp.pi**2), ('1', mp.mpf(1)),
            ('Li2(r)', A(r)), ('Li2(-r)', A(-r)), ('Li2(c)', A(c)), ('Li2(-c)', A(-c)),
            ('Li2((1-c)/2)', A((1-c)/2)), ('Li2((1+c)/2)', A((1+c)/2)),
            ('Li2((1-r)/2)', A((1-r)/2)), ('Li2((1+r)/2)', A((1+r)/2)),
            ('Li2(r/(1+r))', A(r/(1+r))), ('Li2(-r/(1+r))', A(-r/(1+r))),
            ('Li2(2r/(1+r))', A(2*r/(1+r))), ('Li2(1-r)', A(1-r)), ('Li2(1-c)', A(1-c)),
            ('log2^2', L(2)**2), ('logc^2', L(c)**2), ('log(1+c)^2', L(1+c)**2), ('log(1-c)^2', L(1-c)**2),
            ('logr^2', L(r)**2), ('log(1+r)^2', L(1+r)**2), ('log(1-r)^2', L(1-r)**2),
            ('log2logc', L(2)*L(c)), ('log2log(1+c)', L(2)*L(1+c)), ('logclog(1+c)', L(c)*L(1+c)),
            ('logrlog(1+r)', L(r)*L(1+r)), ('logrlog(1-r)', L(r)*L(1-r)), ('logrlog2', L(r)*L(2)),
            ('log(1+r)log(1-r)', L(1+r)*L(1-r)), ('log(1+c)log(1-c)', L(1+c)*L(1-c)), ('log2logr', L(2)*L(r))]

angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6,
          mp.mpf(1)/5, mp.mpf(3)/7, mp.mpf(5)/8, mp.mpf(2)/9, mp.mpf(5)/7, mp.mpf(7)/9,
          mp.mpf(1)/6, mp.mpf(3)/8, mp.mpf(4)/7, mp.mpf(5)/9]
holdout = [mp.mpf(1)/7, mp.mpf(4)/9, mp.mpf(2)/3]
print("evaluating...")
C = [[x for _, x in cand(c)] for c in angles]
names = [n for n, _ in cand(angles[0])]
T = [T_alpha_c(c) for c in angles]
M, N = len(angles), len(names)
# greedy independent column selection (normalized)
A = [mp.matrix([C[i][j] for i in range(M)]) for j in range(N)]
for j in range(N):
    A[j] = A[j]/A[j].norm('2')
sel = []
span = None
for j in range(N):
    v = A[j].copy()
    if span is not None:
        # project out
        for s in sel:
            v = v - (v.T*A[s])[0]*A[s]
    nrm = v.norm('2')
    if nrm > mp.mpf('1e-60'):
        # re-orthogonalize against sel
        sel.append(j)
print(f"independent subset size {len(sel)} of {N}")
# solve T against selected columns (exact M x |sel| least squares via normal eq on independent set)
Asel = mp.matrix(M, len(sel))
for k, j in enumerate(sel):
    for i in range(M):
        Asel[i, k] = C[i][j]
Tv = mp.matrix(T)
G = Asel.T*Asel
b = Asel.T*Tv
q = mp.lu_solve(G, b)
resid = Asel*q - Tv
print("residual:", mp.nstr(resid.norm('inf'), 3))
print("coefficients:")
coef = {}
for k, j in enumerate(sel):
    if abs(q[k]) > mp.mpf('1e-60'):
        coef[names[j]] = q[k]
        print(f"   {names[j]}: {mp.nstr(q[k], 40)}")
# verify holdout
for c in holdout:
    Tc = T_alpha_c(c)
    bmap = {n: x for n, x in cand(c)}
    pred = sum(coef[n]*bmap[n] for n in coef)
    print(f"   holdout c={mp.nstr(c,8)}: diff={mp.nstr(abs(Tc-pred),3)}")
print("STEP20 DONE")
'''
open('step_20.py','w').write(step20)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_20.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
