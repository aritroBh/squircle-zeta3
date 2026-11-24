import mpmath as mp
mp.mp.dps = 60
r2=mp.sqrt(2)
def rho_br(t): return mp.sqrt(1+t*t/(1+r2*t))
def psi(t): return r2*t - mp.log(1+r2*t+t*t)
def H(t):
    rho = rho_br(t)
    return psi(t) - (r2/3)*t**3/(1+r2*t+t*t) - r2*t*(rho-1)**2/rho**2
def Hp3_mp(t):
    den = 1+r2*t+t*t
    q = mp.sqrt((1+r2*t)/den)
    w = 1/(den*(1+q))
    return (r2/3)*(r2+2*t)/den**2 - r2*w*w*t - r2*w*t*(2+r2*t)/(q*den**2)
# H'(t) vs t^3 Hp3(t) via complex-step derivative (exact to machine order)
w = mp.mpf('0')
for k in range(200):
    t = mp.mpf(10)**(mp.mpf(-8) + mp.mpf(9)*k/199)   # t in [1e-8, 10]
    dH = mp.im(H(t+1j*mp.mpf('1e-30')))/mp.mpf('1e-30')
    w = max(w, abs(dH - t**3*Hp3_mp(t))/max(1,abs(dH)))
print(f"R5(b) H'(t) = t^3 Hp3(t) along saddle branch: worst resid (200 pts t in [1e-8,10]) = {mp.nstr(w,6)}")
# H(t) = margin/(N/2) along saddle: margin = lhs - (N t^3/2)(r2/3)/den
mpp=mp.pi; ma=mpp/(2*r2); mbb=mp.e**(1j*mpp/4)
w2 = mp.mpf('0')
for k in range(200):
    t = mp.mpf(10)**(mp.mpf(-8)+mp.mpf(8)*k/199)
    rho = rho_br(t); y = t/rho; N = mpp/(2*y*y); th = mpp/(r2*y); x0=1/(rho*y)
    lhs = mp.re((-1+1j)*ma*x0 + N*mp.log(mbb*x0/(1+mbb*x0))) + th
    margin = lhs - (N*t**3/2)*(r2/3)/(1+r2*t+t*t)
    w2 = max(w2, abs(margin - (N/2)*H(t))/abs(margin))
print(f"      H(t) == margin/(N/2) along saddle: worst rel resid = {mp.nstr(w2,6)}")
# H > 0 for t < 4.5447 ? locate sign change
import mpmath
ts = [mp.mpf(k)/100 for k in range(1, 1001)]
sgn = [H(t) for t in ts]
zc = None
for i in range(len(ts)-1):
    if sgn[i] > 0 and sgn[i+1] <= 0: zc = ts[i]; break
print(f"      H(t) sign change bracket: [{zc}, {zc+mp.mpf('0.01')}]  (D13s claims zero at 4.5447)")
print(f"      H(0.12731) = {mp.nstr(H(mp.mpf('0.12731')),10)} ; H(0.13) = {mp.nstr(H(mp.mpf('0.13')),10)} (certified domain edge, saddle branch)")
