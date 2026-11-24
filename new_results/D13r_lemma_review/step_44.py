import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

def pieces(Nv):
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    y = mp.sqrt(pi/(2*N))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    x0 = 1/(rho*y)
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N); g0 = g(x0)
    om = mp.im(Phi1(x0,N))
    c2 = abs(mp.re(p2))/2
    w = mp.sqrt(mp.log(N)+1.0)
    d = min(w/mp.sqrt(2*c2), x0/3)
    xl, xr = x0-d, x0+d; x2 = x0/2
    sc = N**mp.mpf('-1.5')*mp.e**(-th)
    rt = mp.sqrt(-p2/2)
    tailR = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d + 1j*om/p2))
    tailL = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d - 1j*om/p2))
    TM_exact = abs(g0*mp.e**Phi0*(tailR+tailL))
    def integrand(x): return mp.e**Phi(x,N)*g(x)
    TRt = mp.quad(lambda x: abs(integrand(x)), [xr, x0+8*d, mp.inf])
    TLmt = mp.quad(lambda x: abs(integrand(x)), [x2, xl])
    TLft = mp.quad(lambda x: abs(integrand(x)), [1, x2])
    two_pi = 2/pi
    print(f"N={Nv:g}: (2/pi)*truth/sc: TM={mp.nstr(two_pi*TM_exact/sc,10)} TR={mp.nstr(two_pi*TRt/sc,10)} TLm={mp.nstr(two_pi*TLmt/sc,10)} TLf={mp.nstr(two_pi*TLft/sc,10)}")
    print(f"   buggy certified: TM=0.128969 TR={mp.nstr(1.2533141*0.2361452504,10)} TLm={mp.nstr(1.2533141*1.2907777267,10)}")
pieces(10**10)
pieces(10**12)
# check step_10.py and step_10_J96.py for the same bug
import re
for fn in ['step_10.py','step_10_J96.py','step_10_J192.py']:
    src = open('research/08_theory/D13h_final_cell/'+fn).read()
    for line in src.splitlines():
        if 'ecw' in line and '=' in line and 'def' not in line and 'TMreg' not in line and 'TRreg' not in line and 'TLmreg' not in line:
            print(fn, "::", line.strip())
# corrected final-cell C
corr = 49.56259792396025 + (2/mp.pi)*(mp.pi/2)**mp.mpf('1.5')*(2**mp.mpf('0.25')-1)*(0.10290245010511980+0.23614525040002159+1.29077772666801997)
print("corrected final-cell C (replacing ecw by the right value) =", mp.nstr(corr,12))
