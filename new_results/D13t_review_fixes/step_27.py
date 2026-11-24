import mpmath as mp, random
mp.mp.dps = 60
mpr2=mp.sqrt(2); mpp=mp.pi; ma=mpp/(2*mpr2); mbb=mp.e**(1j*mpp/4)
def Phi(x,N): return (-1+1j)*ma*x + N*mp.log(mbb*x/(1+mbb*x))
def Phi1(x,N): return (-1+1j)*ma + N*(1/x - mbb/(1+mbb*x))

# (1) saddle branch rho(t)=sqrt(1+t^2/(1+rt2 t)) solves Re Phi'(1/t)=0 ?
random.seed(31337)
w1 = mp.mpf('0')
for _ in range(300):
    t = mp.mpf(10)**random.uniform(-14, 2)
    rho = mp.sqrt(1+t*t/(1+mpr2*t)); y = t/rho; N = mpp/(2*y*y)
    w1 = max(w1, abs(mp.re(Phi1(1/t, N))))
print(f"(1) Re Phi'(1/t)=0 on branch rho(t): worst |Re Phi'| = {mp.nstr(w1,6)} over 300 pts t in [1e-14,100]")
# and the branch equation rho^2(1+rt2 t)=1+rt2 t+t^2:
w1b = max(abs((1+t*t/(1+mpr2*t))*(1+mpr2*t)-(1+mpr2*t+t*t)) for t in [mp.mpf(10)**(k/10) for k in range(-140,21)])
print(f"    rho^2(1+rt2 t)=1+rt2 t+t^2 residual: {mp.nstr(w1b,6)}")

# (2) B_true direct vs exact form, dense; also (3) phi identity
def Bform(t):
    rho = mp.sqrt(1+t*t/(1+mpr2*t))
    s1 = -t/((1+mpr2*t)*(1+rho))
    zz = mpr2*t/(2+mpr2*t)
    Ap = 1 - mp.atan(zz)/zz
    A1n = Ap/(t*t)
    B1 = mp.mpf('0.5') + mpr2*s1 + s1*s1 - rho*rho*A1n
    return rho*B1/(2+mpr2*t)
w2 = mp.mpf('0'); w3 = mp.mpf('0'); zzmax = mp.mpf('0')
for k in range(5000):
    t = mp.mpf(10)**(mp.mpf(-14) + mp.mpf(17)*k/4999)   # t in [1e-14, 1e3]
    rho = mp.sqrt(1+t*t/(1+mpr2*t)); y = t/rho; N = mpp/(2*y*y); th = mpp/(mpr2*y); x0 = 1/(rho*y)
    Bt = (mp.im(Phi(x0,N)) - th + mpp/4)/((mpp/mpr2)*y)
    w2 = max(w2, abs(Bt-Bform(t))/abs(Bt))
    zz = mpr2*t/(2+mpr2*t)
    zzmax = max(zzmax, zz)
    phi_id = ma*(1-2*rho)/t + N*mp.atan(zz)      # claimed phi
    phi_true = mp.im(Phi(x0,N)) - th
    w3 = max(w3, abs(phi_id-phi_true))
print(f"(2) B_true exact form vs direct: worst rel resid over 5000 pts t in [1e-14,1e3] = {mp.nstr(w2,6)}")
print(f"(3) zz max over t up to 1e3 = {mp.nstr(zzmax,12)} (<1 always since zz=rt2 t/(2+rt2 t)); phi = a(1-2rho)/t + N atan(zz): worst resid {mp.nstr(w3,6)}")
