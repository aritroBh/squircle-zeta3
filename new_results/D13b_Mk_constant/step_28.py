code = r'''
import mpmath as mp
mp.mp.dps = 60

def B(z):
    return 1 - mp.exp(z)/(1+mp.sqrt(1+mp.exp(2*z)))

def Discf(s, r):
    zs = s*1j*mp.pi/2
    dlt = mp.mpf('1e-42')
    return B(zs + r + 1j*dlt) - B(zs + r - 1j*dlt)

def Jhat(s, u, N, R):
    # (1/2pii) ∫_0^R Disc(r) (1+(zs+r)/u)^{-N} dr
    zs = s*1j*mp.pi/2
    f = lambda r: Discf(s, r)*(1+(zs+r)/u)**(-N)
    return mp.quad(f, [0, R])/(2*mp.pi*1j)

def saddle(N, s):
    zs = s*1j*mp.pi/2
    us = (-zs + mp.sqrt(zs**2 + 4*N*zs))/2
    # phi''(us) exactly:  phi'' = -N zs (2u+zs)/(u(us+zs))^2, and u(us+zs)=N zs
    phipp = -(2*us+zs)/(N*zs)
    return us, phipp

def Iu_exact(N, s, npts=60):
    us, phipp = saddle(N, s)
    sig = 1/mp.sqrt(-phipp)         # complex; use abs
    T = 8*abs(sig)
    # steepest direction: dir with phipp*dir^2 real negative
    argd = (mp.pi - mp.arg(phipp))/2
    dirv = mp.e**(1j*argd)
    R = abs(us)*(mp.e**(45*mp.log(10)/N))   # where (1+R/|us|)^{-N} ~ 1e-45
    R = min(R, mp.mpf('2.0'))
    def integrand(t):
        u = us + t*dirv
        return mp.e**(-u)*Jhat(s, u, N, R)/u*dirv
    pts = [T*mp.cos(mp.pi*i/npts) for i in range(npts+1)]  # Chebyshev
    # use tanh-sinh via mp.quad on [-T,T]
    return mp.quad(integrand, [-T, 0, T])

def Iu_analytic(N, s, beta_s):
    zs = s*1j*mp.pi/2
    us, phipp = saddle(N, s)
    phi_s = -us - N*mp.log(1+zs/us)
    pre = (2*s*mp.sqrt(2)/(2*mp.pi*1j))*mp.gamma(mp.mpf('1.5'))*N**mp.mpf('-1.5')*mp.sqrt(mp.pi)*mp.sqrt(zs)*mp.sqrt(N)
    return pre*mp.e**phi_s*(1 + beta_s/mp.sqrt(N))

for N in [mp.mpf(401), mp.mpf(1601)]:
    zs = 1j*mp.pi/2
    beta_s = 1/(2*mp.sqrt(zs)) + mp.sqrt(zs)/2 - zs**mp.mpf('1.5')/12
    Ie = Iu_exact(N, 1)
    Ia0 = Iu_analytic(N, 1, 0)
    Ia1 = Iu_analytic(N, 1, beta_s)
    r = Ie/Ia0 - 1
    print(f"N={int(N)}: ratio_exact/leading-1 = {mp.nstr(r,12)}")
    print(f"   => implied beta_s = sqrt(N)*r = {mp.nstr(mp.sqrt(N)*r,12)}  (up to O(1/sqrt N))")
    print(f"   analytic beta_s = {mp.nstr(beta_s,12)}")
    print(f"   ratio exact/with-beta = {mp.nstr(Ie/Ia1,10)}")
'''
open('step_16.py','w').write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_16.py"], capture_output=True, text=True, timeout=1800)
print(r.stdout); print(r.stderr[-1200:] if r.stderr else "")
