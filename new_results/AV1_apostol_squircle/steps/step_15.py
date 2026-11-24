
# step_15.py : TASK 4 special points: verify proposed closed forms at 90 digits.
import mpmath as mp
mp.mp.dps = 100

def A_series(p, s=0, tol=None):
    p = mp.mpf(p); s = mp.mpf(s); tol = tol or mp.mpf('1e-95')
    tot = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma((m+s)/p)**2 / (p**2 * mp.gamma(1 + 2*(m+s)/p))
        tot += t
        if m > 30 and abs(t) < tol: break
    return tot

def B_series(p, tol=None):
    p = mp.mpf(p); tol = tol or mp.mpf('1e-95')
    tot = mp.mpf(0); n = -1
    while True:
        n += 1
        t = mp.gamma((2*n+1)/p)**2 / (p**2 * mp.gamma(1 + 2*(2*n+1)/p))
        tot += t
        if n > 30 and abs(t) < tol: break
    return tot

z2 = mp.zeta(2); phi = (1+mp.sqrt(5))/2
print('=== p=2 quarter disc ===')
A2 = A_series(2)
print('A(2) =', mp.nstr(A2, 30), ' 7*zeta2/12 =', mp.nstr(7*z2/12, 30), ' err', mp.nstr(abs(A2-7*z2/12),3))
# binomial ingredients
S1 = mp.nsum(lambda k: mp.binomial(2*k,k)/(16**k*(2*k-1)), [1,mp.inf])
S2 = mp.nsum(lambda k: mp.binomial(2*k,k)/(16**k*(2*k-1)**2), [1,mp.inf])
print('sum C(2k,k)/(16^k(2k-1)) =', mp.nstr(S1,25), ' (2-sqrt3)/2 =', mp.nstr((2-mp.sqrt(3))/2,25), ' err', mp.nstr(abs(S1-(2-mp.sqrt(3))/2),3))
print('sum C(2k,k)/(16^k(2k-1)^2) =', mp.nstr(S2,25), ' pi/12-(2-sqrt3)/2 =', mp.nstr(mp.pi/12-(2-mp.sqrt(3))/2,25), ' err', mp.nstr(abs(S2-(mp.pi/12-(2-mp.sqrt(3))/2)),3))

print('=== p=1/2 ===')
Ah = A_series(mp.mpf('0.5'))
cand = 2*z2/3 - 4*mp.log(phi)**2
print('A(1/2) =', mp.nstr(Ah,30), ' 2z2/3 - 4 log^2 phi =', mp.nstr(cand,30), ' err', mp.nstr(abs(Ah-cand),3))
# C = iint_tri dxdy/(1+xy) = 2 log^2 phi ?
C = mp.quad(lambda x: mp.log(1+x*(1-x))/x, [0, mp.mpf('0.5'), 1])
print('C =', mp.nstr(C,25), ' 2 log^2 phi =', mp.nstr(2*mp.log(phi)**2,25), ' err', mp.nstr(abs(C-2*mp.log(phi)**2),3))

print('=== B(2), B(1), B(1/2) ===')
B2 = B_series(2); print('B(2) =', mp.nstr(B2,30), ' pi^2/12 =', mp.nstr(mp.pi**2/12,30), ' err', mp.nstr(abs(B2-mp.pi**2/12),3))
B1 = B_series(1); print('B(1) =', mp.nstr(B1,30), ' z2/6+log^2 phi =', mp.nstr(z2/6+mp.log(phi)**2,30), ' err', mp.nstr(abs(B1-z2/6-mp.log(phi)**2),3))
Bh = B_series(mp.mpf('0.5')); print('B(1/2) =', mp.nstr(Bh,30))
print('B(1/2) PSLQ {1,pi,log2,log3,z2,log phi}:', mp.pslq([Bh, mp.mpf(1), mp.pi, mp.log(2), mp.log(3), mp.zeta(2), mp.log(phi)], tol=mp.mpf('1e-40'), maxcoeff=10**7, maxsteps=10000))
print('B(1/2) PSLQ {z2, log^2 phi, log phi, log2^2}:', mp.pslq([Bh, mp.zeta(2), mp.log(phi)**2, mp.log(phi), mp.log(2)**2], tol=mp.mpf('1e-40'), maxcoeff=10**7, maxsteps=10000))

print('=== p=3/2 ===')
A32 = A_series(mp.mpf('1.5'))
print('A(3/2) =', mp.nstr(A32,30))
print('A(3/2) PSLQ {1,pi,log2,z2,z3,pi/rt3,G(1/3)^3/sqrt3}:', mp.pslq([A32, mp.mpf(1), mp.pi, mp.log(2), mp.zeta(2), mp.zeta(3), mp.pi/mp.sqrt(3), mp.gamma(mp.mpf(1)/3)**3/mp.sqrt(3)], tol=mp.mpf('1e-38'), maxcoeff=10**7, maxsteps=10000))

print('=== A_1(1) = iint_tri xy/(1-xy) ===')
A11 = A_series(1, s=1)
print('A_1(1) =', mp.nstr(A11,30))
print('PSLQ {1, pi/rt3, z2, log2}:', mp.pslq([A11, mp.mpf(1), mp.pi/mp.sqrt(3), mp.zeta(2), mp.log(2)], tol=mp.mpf('1e-40'), maxcoeff=10**7, maxsteps=10000))
print('PSLQ {1, z2, pi}:', mp.pslq([A11, mp.mpf(1), mp.zeta(2), mp.pi], tol=mp.mpf('1e-40'), maxcoeff=10**7, maxsteps=10000))
