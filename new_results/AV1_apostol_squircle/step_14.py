code = r'''
# step_11.py : TRIANGLE END. A(1)=zeta(2)/3 check; expansion of A in tau = 1-1/p at p=1.
# A(1/(1-tau)) = sum_k c_k tau^k,  c_k = (-1)^k/k! sum_{m>=1} m^{k-2} f^{(k)}(m),
# f(z)=Gamma(1+z)^2/Gamma(1+2z), log-derivatives g^{(j)} via polygamma; Bell polynomials.
import mpmath as mp
mp.mp.dps = 100

def f(m):
    m = int(m)
    return mp.gamma(m+1)**2/mp.gamma(2*m+1)

def gder(m, j):  # g^{(j)}(m) = 2 psi^{(j-1)}(1+m) - 2^j psi^{(j-1)}(1+2m)
    return 2*mp.psi(j-1, m+1) - mp.mpf(2)**j*mp.psi(j-1, 2*m+1)

from functools import lru_cache
def bell(fderivs):
    # given [g1..gk], return [f1..fk] complete Bell polynomial values f^{(k)}/f
    k = len(fderivs)
    B = [[mp.mpf(0)]*(k+1) for _ in range(k+1)]
    B[0][0] = mp.mpf(1)
    for n in range(1, k+1):
        for j in range(1, n+1):
            B[n][j] = sum(mp.binomial(n-1,i-1)*fderivs[i-1]*B[n-i][j-1] for i in range(1, n-j+2))
    return [B[n][1] for n in range(1,k+1)]

def fder(m, k):
    gs = [gder(m, j) for j in range(1, k+1)]
    Bs = bell(gs)
    fm = f(m)
    return [fm*b for b in Bs]

K = 6
cs = [mp.mpf(0)]*(K+1)
M = 500
for m in range(1, M+1):
    fm = f(m)
    if fm < mp.mpf('1e-110') and m > 50: break
    fds = fder(m, K)
    for k in range(0, K+1):
        term = (mp.mpf(m)**(k-2)) * (fds[k-1] if k>0 else fm)
        cs[k] += (-1)**k/mp.factorial(k)*term
for k in range(K+1):
    print(f'c_{k} =', mp.nstr(cs[k], 40))
print('zeta2/3 =', mp.nstr(mp.zeta(2)/3, 40))
print('|c0 - zeta2/3| =', mp.nstr(abs(cs[0]-mp.zeta(2)/3),3))

# PSLQ pools
rt3 = mp.sqrt(3)
def report(name, val, basis, names, maxcoeff=10**6):
    rel = mp.pslq([val]+basis, tol=mp.mpf('1e-40'), maxcoeff=maxcoeff, maxsteps=10000)
    print(f'{name}: PSLQ {rel} vs {names}')
    return rel
print('--- c1 pools')
report('c1', cs[1], [mp.pi/rt3, mp.log(2), mp.mpf(1)], ['pi/sqrt3','log2','1'])
report('c1', cs[1], [mp.pi/rt3, mp.log(2), mp.zeta(2), mp.mpf(1)], ['pi/sqrt3','log2','z2','1'])
print('--- c2 pools')
report('c2', cs[2], [mp.zeta(2), mp.pi**2/rt3, mp.zeta(3), mp.pi/rt3, mp.log(2)**2, mp.mpf(1)], ['z2','pi^2/sqrt3','z3','pi/sqrt3','log2^2','1'])
print('--- c3 pools')
report('c3', cs[3], [mp.zeta(3), mp.pi*mp.zeta(2)/rt3, mp.pi**3/rt3, mp.zeta(2), mp.pi/rt3, mp.log(2)*mp.zeta(2), mp.mpf(1)], ['z3','pi*z2/sqrt3','pi^3/sqrt3','z2','pi/sqrt3','log2*z2','1'])
'''
open('steps/step_11.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_11.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
open('progress.log','a').write('step_11: triangle-end tau-expansion coefficients c0..c6 computed, PSLQ\n')
