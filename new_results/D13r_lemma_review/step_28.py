import sys
sys.path.insert(0, 'research/08_theory/D13g_all_N')
from flint import arb, acb, ctx
ctx.dps = 60
import mpmath as mp
mp.mp.dps = 50
pi = arb.pi(); rt2 = arb(2).sqrt(); a = pi/(2*rt2); sqpi = pi.sqrt()
b = (acb(0,1)*pi/4).exp(); I1 = acb(0,1)
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
def Hc(s): return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
def gF(x): return Hc(b*x)/x
def h3encl(T):
    T = arb(T).upper()
    return iv(rt2/3/(1+rt2*T+T*T), rt2/3)
def G1bound(x):
    t = x/(2*rt2)
    return arb("0.75")*(1/t**2 + 1/t) + 1/x

def zbox(y1, y2, nsub=4):
    """replicate step_38's z computation over the (y,sig) box, return arb complex box for z"""
    y = iv(y1, y2)
    zb = None
    for k in range(nsub):
        sig = iv(arb(k)/nsub, arb(k+1)/nsub)
        rho = 1 + y*y*sig
        y0 = rho*y; x0 = 1/y0
        N = pi/(2*y*y)
        x0lo = arb(x0.lower()); x0hi = arb(x0.upper())
        den1 = 1+rt2*rho*y+rho*rho*y*y
        Omn = rt2*rho + rho*rho*y - 2*y*sig - y**3*sig**2
        om = (pi/(2*rt2))*y*Omn/den1
        mS = (N*y0**3/2)*h3encl(y0)
        xref = arb(x0.mid()); hx = (x0hi-x0lo)/2
        rg = arb(G1bound(x0lo).upper())*hx
        g0box = gF(xref) * acb(arb(0,rg), arb(0,rg)).exp()
        t = rho*y
        numR = t**3 + 2*rt2*t**2 + 4*t + rt2
        D2 = t**2 + 2*rt2*t + 4
        invp2 = -(2/(pi*rho**3*y))*acb(numR, t+rt2)/D2
        zz = rt2*rho*y/(2+rt2*rho*y); z2 = zz*zz
        A1m = iv(0, (z2/(1+z2)).upper())
        T2n = 2*y*y*sig - rt2*rho*y - 2*rho*A1m*y
        Bh = -y*sig/(2*rho) + T2n/(2*y*(2+rt2*rho*y))
        phireg = (pi/rt2)*Bh
        eiphi = (I1*phireg).exp()
        z = (2*N*I1/pi**arb("1.5"))*g0box*mS.exp()*((-2*pi*invp2).sqrt())*((om**2*invp2/2).exp())*eiphi
        zb = z if zb is None else zb.union(z) if hasattr(zb,'union') else zb
        # acb has no union in older flint; collect manually
        if k==0: zs=[z]
        else: zs.append(z)
    return zs

mpp = mp.pi; msqrt2 = mp.sqrt(2); ma = mpp/(2*msqrt2); mbb = mp.e**(1j*mpp/4)
def mHc(s): return mp.sqrt(mpp)/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def mPhi(x, N): return (-1+1j)*ma*x + N*mp.log(mbb*x/(1+mbb*x))
def mPhi2(x, N): return -N*(1+2*mbb*x)/(x**2*(1+mbb*x)**2)
def z_true(N):
    y = mp.sqrt(mpp/(2*N)); th = mpp/(msqrt2*y)
    f = lambda r: r**2*(1+msqrt2*r*y)/(1+msqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-40'))
    x0 = 1/(rho*y)
    Phi0 = mPhi(x0,N); p2 = mPhi2(x0,N)
    om = mp.im((-1+1j)*ma + N*(1/x0 - mbb/(1+mbb*x0)))
    MV = (mHc(mbb*x0)/x0)*mp.e**Phi0*mp.sqrt(-2*mpp/p2)*mp.e**(om**2/(2*p2))
    Lt = -1j*(mpp**mp.mpf('1.5')/(2*N))*mp.e**(-th)*mp.e**(1j*th)
    return MV/Lt

def acb_contains(zbox_list, zt):
    """does any of the acb boxes contain zt (as a small mpc box)?"""
    for z in zbox_list:
        re_ok = z.real.lower() <= arb(float(mp.re(zt))) <= z.real.upper()
        im_ok = z.imag.lower() <= arb(float(mp.im(zt))) <= z.imag.upper()
        if re_ok and im_ok: return True
    return False

for (N1,N2) in [(100,150),(150,225),(225,340),(340,400)]:
    y1 = mp.sqrt(mpp/(2*N2)); y2 = mp.sqrt(mpp/(2*N1))   # y decreases with N
    zs = zbox(float(y1)*(1-1e-9), float(y2)*(1+1e-9), nsub=4)
    miss = 0; tot = 0; marg = mp.mpf('-9')
    for k in range(41):
        Nv = N1 + (N2-N1)*k/40
        zt = z_true(mp.mpf(Nv))
        tot += 1
        ok = False
        for z in zs:
            # distance from zt to box
            dx = max(float(z.real.lower()) - float(mp.re(zt)), float(mp.re(zt)) - float(z.real.upper()), 0.0)
            dy = max(float(z.imag.lower()) - float(mp.im(zt)), float(mp.im(zt)) - float(z.imag.upper()), 0.0)
            d = (dx**2+dy**2)**0.5
            if d == 0.0: ok=True
            marg = max(marg, -mp.mpf(d))
        if not ok:
            miss += 1
    print(f"cell [{N1},{N2}]: z_true OUTSIDE step_38 z-box for {miss}/{tot} sample N")
