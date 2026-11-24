# step_07.py -- TASK 3: growth analysis, radius, singularity asymptotics, optimal truncation
import mpmath as mp, pickle
mp.mp.dps = 80
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']; INT=pickle.load(open('/tmp/Nvals.pkl','rb'))['INT']
pi=mp.pi
print("k : N_k/prev  |N_k|^{1/k}  k/(e*pi)  sqrt(|N_{k+2}/N_k|)*pi/(k+1)  odd-norm N_k*pi^{k-1}/(2(k-2)!)")
for k in range(10,41):
    rt=abs(N[k])**(1/k); q=mp.sqrt(abs(N[k+2]/N[k]))*pi/(k+1)
    norm = N[k]*pi**(k-1)/(2*mp.factorial(k-2)) if k%2==1 else mp.nan
    print(f"{k:2d} {mp.nstr(N[k+1]/N[k],8):>10} {mp.nstr(rt,8)} {mp.nstr(k/(mp.e*pi),8)} {mp.nstr(q,8)}  {mp.nstr(norm,8) if k%2==1 else ''}")
# integrand singularity asymptotics: [h^m]psi2 ~ -(2/s) Re(1+i s/pi)^{m+1}
def psi2(s,h): return (1+mp.e**(-s*h/(1-h)))**(-1-h)
def coeff_psi2(s,m):
    rho=0.8*pi/mp.sqrt(s**2+pi**2)
    return mp.quad(lambda t: psi2(s,rho*mp.e**(1j*t))/(rho*mp.e**(1j*t))**m,[0,2*pi])/(2*pi)
for s in [10.0,30.0]:
    for m in [9,10]:
        cm=coeff_psi2(s,m); asym=-(1/s)*2*mp.re((1+1j*s/pi)**(m+1))
        print(f"s={s} m={m}: [h^m]psi2={mp.nstr(cm,10)} asym={mp.nstr(asym,10)} ratio={mp.nstr(cm/asym,6)}")
# optimal truncation at p=10
mp.mp.dps=110
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']
def pi_cosh(p):
    p=mp.mpf(p); h=1/p
    f=lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f,[0,mp.inf])
true=pi_cosh(10)/4; h=mp.mpf('0.1'); s=mp.mpf(1); best=(0,mp.mpf(1))
for k in range(0,42):
    s+=N[k]*h**(k+1); err=abs(s-true)
    if err<best[1]: best=(k,err)
print("optimal truncation p=10: best k=",best[0]," err=",mp.nstr(best[1],3),
      " pi/h=",mp.nstr(pi/h,5)," e^(-pi/h)=",mp.nstr(mp.e**(-pi/h),3))
