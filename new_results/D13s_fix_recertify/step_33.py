
import os, subprocess
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'

zcode = r'''
import sys, math
import mpmath as mp
mp.mp.dps = 50
sys.path.insert(0, '.')
from flint import arb, acb, ctx
ctx.dps = 60
from step_38 import iv, pi, rt2, gF, G1bound, h3encl
I1 = acb(0,1); b = (I1*pi/4).exp()

def zbox_fixed(y1, y2, nsub=4):
    """replicate the FIXED step_38 z computation over the (y,sig) box"""
    y = iv(y1, y2); zs = []
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
        zz = rt2*rho*y/(2+rt2*rho*y)
        # FIXED A1m
        zz_lo = arb(zz.lower()); zz_up = arb(zz.upper())
        fA1_lo = arb((1 - zz_lo.atan()/zz_lo).lower())
        fA1_up = arb((1 - zz_up.atan()/zz_up).upper())
        A1m = iv(fA1_lo/arb(y.upper()), fA1_up/arb(y.lower()))
        T2n = 2*y*y*sig - rt2*rho*y - 2*rho*A1m*y
        Bh = -y*sig/(2*rho) + T2n/(2*y*(2+rt2*rho*y))
        phireg = (pi/rt2)*Bh
        eiphi = (I1*phireg).exp()
        z = (2*N*I1/pi**arb("1.5"))*g0box*mS.exp()*((-2*pi*invp2).sqrt())*((om**2*invp2/2).exp())*eiphi
        zs.append(z)
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

total_miss = 0; total = 0; worst_margin = mp.mpf('9')
for (N1,N2) in [(100,150),(150,225),(225,340),(340,400),(400,600),(600,900),(1e4,1.5e4),(1e6,1.5e6),(1e9,1.5e9)]:
    y1 = math.sqrt(math.pi/(2*N2)); y2 = math.sqrt(math.pi/(2*N1))
    zs = zbox_fixed(y1*(1-1e-9), y2*(1+1e-9), nsub=4)
    miss = 0; marg = mp.mpf('9')
    for k in range(41):
        Nv = N1 + (N2-N1)*k/40
        zt = z_true(mp.mpf(Nv))
        total += 1
        ok = False
        bestd = mp.mpf('9')
        for z in zs:
            dx = max(float(z.real.lower()) - float(mp.re(zt)), float(mp.re(zt)) - float(z.real.upper()), 0.0)
            dy = max(float(z.imag.lower()) - float(mp.im(zt)), float(mp.im(zt)) - float(z.imag.upper()), 0.0)
            d = (dx**2+dy**2)**0.5
            bestd = min(bestd, mp.mpf(d))
            if d == 0.0: ok=True
        if not ok:
            miss += 1; total_miss += 1
        marg = min(marg, -bestd)
    print(f"cell [{N1},{N2}]: misses {miss}/41   (min signed margin {mp.nstr(marg,6)})")
print(f"TOTAL: {total_miss} misses out of {total} sampled N")
'''
open(os.path.join(FIX,'verify_zbox.py'),'w').write(zcode)
r = subprocess.run([PY, 'verify_zbox.py'], cwd=FIX, capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-1200:] if r.stderr else "")
