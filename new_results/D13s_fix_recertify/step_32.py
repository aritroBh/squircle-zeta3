
import os, subprocess
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'

# (a) Bh exactness with A1 = A1'/y (review step_29 claim, independent check)
code_a = r'''
import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
random.seed(43)
w = mp.mpf('0')
for _ in range(400):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    sig = (rho-1)/y**2; x0 = 1/(rho*y); t = rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t)
    Ap = 1 - mp.atan(zz)/zz
    A1 = Ap/y                      # the value that makes Bh exact
    T2n = 2*y*y*sig - sqrt2*rho*y - 2*rho*A1*y
    Bh = -y*sig/(2*rho) + T2n/(2*y*(2+sqrt2*t))
    w = max(w, abs(phi_true - (pi/sqrt2)*Bh))
print("phi = (pi/rt2) Bh[A1'/y]: worst residual over 400 pts =", mp.nstr(w,10))
'''
open(os.path.join(FIX,'verify_Bh.py'),'w').write(code_a)
r = subprocess.run([PY, 'verify_Bh.py'], cwd=FIX, capture_output=True, text=True, timeout=1200)
print(r.stdout); print(r.stderr[-1200:] if r.stderr else "")
