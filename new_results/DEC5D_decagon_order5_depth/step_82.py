
import mpmath as mp
mp.mp.dps = 250
z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
ZV = [mp.mpf(1), z2, z3, z2**2, z5, z2*z3]
ZN = ['1','z2','z3','z2^2','z5','z2z3']
ZW = [0,2,3,4,5,5]
v = -mp.mpf(5)/1944
cand = [k for k in range(6) if ZW[k] <= 5]
vec = [v] + [ZV[k] for k in cand]
rel = mp.pslq(vec, tol=mp.mpf('1e-200'), maxcoeff=10**25, maxsteps=10000)
print("(5,5,0) rel:", rel)
# try with smaller tol
rel2 = mp.pslq(vec, tol=mp.mpf('1e-150'), maxcoeff=10**25, maxsteps=10000)
print("(5,5,0) rel2:", rel2)
# a weight-5 mixed one: (0,2,3): int (log t)^2 (log(1-t))^3
f = lambda x: mp.log(x)**2*mp.log(1-x)**3
v2 = mp.quad(f, [0,1], method='tanh-sinh')
vec2 = [v2] + ZV
rel3 = mp.pslq(vec2, tol=mp.mpf('1e-200'), maxcoeff=10**25, maxsteps=10000)
print("(0,2,3) rel3:", rel3)
