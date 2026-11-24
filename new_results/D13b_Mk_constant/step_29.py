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
    zs = s*1j*mp.pi/2
    f = lambda r: Discf(s, r)*(1+(zs+r)/u)**(-N)
    return mp.quad(f, [0, R])/(2*mp.pi*1j)

def run(N, s, h=mp.mpf('0.18'), T=mp.mpf('7')):
    zs = s*1j*mp.pi/2
    us = (-zs + mp.sqrt(zs**2 + 4*N*zs))/2
    phi  = lambda u: -u - N*mp.log(1+zs/u)
    phip = lambda u: -1 + N*zs/(u*(u+zs))
    phis = phi(us)
    phipp = -(2*us+zs)/(N*zs)
    phippp = lambda u: -N*zs*( -1/(u**2*(u+zs)) - 1/(u*(u+zs)**2) ) # d/du[N zs/(u(u+zs))]
    # series of u(t) at 0: u = us + a1 t + a2 t^2/2 + a3 t^3/6, phi = phis - t^2
    p2 = phipp
    p3 = -N*zs*( -1/(us**2*(us+zs)) - 1/(us*(us+zs)**2) )
    p4 = -N*zs*( 2/(us**3*(us+zs)) + 1/(us**2*(us+zs)**2) + 1/(us**2*(us+zs)**2) + 2/(us*(us+zs)**3) )
    a1 = mp.sqrt(-2/p2)
    # choose branch with Re(a1)>0-ish (direction of increasing t)  -- pick arg consistent w/ steepest
    a2 = -p3*a1**3/(3*p2*a1)  # from p2 a1 a2 + p3 a1^3/6 = 0
    a2 = -p3*a1**2/(3*p2)
    a3 = -(p4*a1**4 + 3*p2*a2**2 + 3*p3*a1**2*a2)/(3*p2*a1)/2
    R = min(abs(us)*(mp.e**(45*mp.log(10)/N)), mp.mpf('2.5'))
    def u_of_t(t, uprev):
        if abs(t) < mp.mpf('0.4'):
            return us + a1*t + a2*t**2/2 + a3*t**3/6
        u = uprev
        for _ in range(40):
            F = phi(u) - (phis - t**2)
            u = u - F/phip(u)
        return u
    def dudt(t, u):
        if abs(t) < mp.mpf('0.4'):
            return a1 + a2*t + a3*t**2/2
        return -2*t/phip(u)
    tot = mp.mpc(0)
    n = int(T/h)
    # integrate t in [-T, T]: two sweeps from 0 outward for Newton continuation
    for sgn in [+1, -1]:
        uprev = u_of_t(sgn*mp.mpf('0.4'), None) if False else None
        ts = [sgn*h*i for i in range(0, n+1)]
        us_path = []
        u0 = us + a1*ts[0] + a2*ts[0]**2/2 + a3*ts[0]**3/6
        us_path.append(u0)
        for i in range(1, len(ts)):
            t = ts[i]
            if abs(t) < mp.mpf('0.4'):
                us_path.append(us + a1*t + a2*t**2/2 + a3*t**3/6)
            else:
                us_path.append(u_of_t(t, us_path[-1]))
        ssum = mp.mpc(0)
        for i, t in enumerate(ts):
            u = us_path[i]
            w = mp.e**(phis - t**2 - u)*Jhat(s, u, N, R)/u*dudt(t, u)
            wt = w*h*(mp.mpf('0.5') if (i==0 or i==len(ts)-1) else 1)
            ssum += wt
        tot += ssum
    return tot

def Iu_leading(N, s):
    zs = s*1j*mp.pi/2
    us = (-zs + mp.sqrt(zs**2 + 4*N*zs))/2
    phi_s = -us - N*mp.log(1+zs/us)
    pre = (2*s*mp.sqrt(2)/(2*mp.pi*1j))*mp.gamma(mp.mpf('1.5'))*N**mp.mpf('-1.5')*mp.sqrt(mp.pi)*mp.sqrt(zs)*mp.sqrt(N)
    return pre*mp.e**phi_s

for N in [mp.mpf(401), mp.mpf(1601)]:
    Ie = run(N, 1)
    Il = Iu_leading(N, 1)
    r = Ie/Il - 1
    print(f"N={int(N)}: exact/leading - 1 = {mp.nstr(r,12)}   implied beta_s = {mp.nstr(mp.sqrt(N)*r,12)}")
zs = 1j*mp.pi/2
print("analytic beta_s =", mp.nstr(1/(2*mp.sqrt(zs)) + mp.sqrt(zs)/2 - zs**mp.mpf('1.5')/12, 14))
'''
open('step_17.py','w').write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_17.py"], capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-1200:] if r.stderr else "")
