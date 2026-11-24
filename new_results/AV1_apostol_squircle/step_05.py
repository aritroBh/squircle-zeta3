code = r'''
# step_04.py : square-end coefficients to 40 digits, three ways; PSLQ.
# d1 = int log(-Phi/u); d2 = (1/2)int L; d3 = (1/24)int(Phi^2-u^2); d5 = -(1/2880)int(Phi^4-u^4);
# d7 = (1/181440)int(Phi^6-u^6).   L(u)=log(1-e^-u), Phi=-u+L.
import mpmath as mp
mp.mp.dps = 105

def A_series(p, tol=mp.mpf('1e-101')):
    p = mp.mpf(p); s = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma(1 + m/p)**2 / (m**2 * mp.gamma(1 + 2*m/p))
        s += t
        if m > 20 and abs(t) < tol: break
    return s

def L(u): return mp.log(1 - mp.exp(-u))
def Phi(u): return -u + L(u)
pts = [mp.mpf('1e-12'), mp.mpf('1e-8'), mp.mpf('1e-5'), mp.mpf('1e-3'), mp.mpf('0.05'), mp.mpf('0.3'), 1, 3, 8, 20, mp.inf]

d1 = mp.quad(lambda u: mp.log(-Phi(u)/u), pts)
d2 = mp.quad(lambda u: L(u)/2, pts)
d3 = mp.quad(lambda u: (Phi(u)**2 - u**2)/24, pts)
d5 = -mp.quad(lambda u: (Phi(u)**4 - u**4), pts)/2880
d7 = mp.quad(lambda u: (Phi(u)**6 - u**6), pts)/181440
print('d1 =', mp.nstr(d1,45))
print('d2 =', mp.nstr(d2,45), '  -zeta2/2 =', mp.nstr(-mp.zeta(2)/2,45), ' err', mp.nstr(abs(d2+mp.zeta(2)/2),3))
print('d3 =', mp.nstr(d3,45), '  zeta3/6  =', mp.nstr(mp.zeta(3)/6,45), ' err', mp.nstr(abs(d3-mp.zeta(3)/6),3))
print('d5 =', mp.nstr(d5,45))
print('d7 =', mp.nstr(d7,45))

# cross-check via Richardson on the series: D(eps)=zeta2-A; fit D/eps = d1+d2*eps+d3*eps^2+0*eps^3+d5*eps^4+...
ps = list(range(40, 401, 8))
xs = [mp.mpf(1)/p for p in ps]
ys = [(mp.zeta(2)-A_series(p))*p for p in ps]
deg = 9
# least squares
n = len(xs); M = mp.matrix(n, deg+1)
for i,x in enumerate(xs):
    for j in range(deg+1): M[i,j] = x**j
b = mp.matrix(n,1)
for i,y in enumerate(ys): b[i,0] = y
# normal equations
Mt = M.T
sol = mp.lu_solve(Mt*M, Mt*b)
print('series-fit d1..d10:', [mp.nstr(sol[j,0],12) for j in range(deg+1)])

# PSLQ d5 against weight-5 basis {zeta5, zeta2*zeta3}
for name, val, basis in [('d5', d5, [mp.zeta(5), mp.zeta(2)*mp.zeta(3)]),
                          ('d1', d1, None)]:
    if basis:
        rel = mp.pslq([val]+basis, tol=mp.mpf('1e-38'), maxcoeff=mp.mpf('1e12'))
        print(f'PSLQ {name} vs weight-5 basis:', rel)
# PSLQ d7 against weight-8 basis {zeta8, zeta3*zeta5, zeta2*zeta3^2, zeta(3,5)}
z35 = mp.zeta(3,5) if hasattr(mp,'zeta') else None
try:
    z35 = mp.zeta(mp.mpf(3),mp.mpf(5))
except Exception as e:
    print('no mzeta in this mpmath:', e)
    z35 = None
print('zeta8 =', mp.nstr(mp.zeta(8),20), 'z3z5 =', mp.nstr(mp.zeta(3)*mp.zeta(5),20), 'z2z3^2 =', mp.nstr(mp.zeta(2)*mp.zeta(3)**2,20), 'z(3,5)=', mp.nstr(z35,20) if z35 else None)
'''
open('steps/step_04.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_04.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')
open('progress.log','a').write('step_04: d1..d7 at 45 digits via corner integrals; PSLQ attempts\n')
