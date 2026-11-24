import sys, math
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')
iv, h3encl, pi, rt2 = st.iv, st.h3encl, st.pi, st.rt2
sq2 = rt2; a2 = pi*sq2/4
import mpmath as mp
mp.mp.dps = 60; r2 = mp.sqrt(2); mpp = mp.pi

def truth_pt(y, sig):
    y = arb(y); sig = arb(sig); rho = 1+y*y*sig; t = rho*y; N = pi/(2*y*y)
    w = 1 - (1+sq2*t+t*t)*(-sq2*t).exp(); psi = -(1-w).log()
    return (N/2)*psi - (a2/t)*(rho-1)**2
def coded_mS(y1,y2,s1,s2):
    y = iv(y1,y2); sig = iv(s1,s2); rho = 1+y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)
def saddle_sig(yv):
    yv=mp.mpf(yv); f=lambda tt: tt/mp.sqrt(1+tt**2/(1+r2*tt))-yv
    t=mp.findroot(f, yv*(1+yv**2/2)); rho=mp.sqrt(1+t**2/(1+r2*tt if False else 1+r2*t))
    return (rho-1)/yv**2

cells = []; N0=100.0
while N0 < 400*(1-1e-12):
    N1=min(N0*1.05,400.0); cells.append((N0,N1)); N0=N1
w_lo=None; w_hi=None
for (n1,n2) in cells:
    y1=math.sqrt(math.pi/(2*n2)); y2=math.sqrt(math.pi/(2*n1))
    mS=coded_mS(y1,y2,0.25,0.5); mlo=arb(mS.lower()); mup=arb(mS.upper())
    for yv in [y1,(y1+y2)/2,y2]:
        sg=float(saddle_sig(yv)); tv=truth_pt(yv,sg)
        m_lo=arb(tv.lower())-mlo; m_hi=mup-arb(tv.upper())
        if w_lo is None or arb(m_lo.upper())<arb(w_lo.upper()): w_lo=m_lo
        if w_hi is None or arb(m_hi.upper())<arb(w_hi.upper()): w_hi=m_hi
print(f"saddle worst lower margin = {arb(w_lo.lower())}")
print(f"saddle worst upper margin = {arb(w_hi.lower())}")

# sign structure of G(y,1) on (0, y100]
def psi(t): return r2*t - mp.log(1+r2*t+t*t)
def G(y, sg):
    rho=1+y*y*sg; t=rho*y
    return psi(t)-(r2/3)*t**3/(1+r2*t+t*t)-r2*t*(rho-1)**2/rho**2
y100 = mp.sqrt(mpp/200)
prev=None; crossings=[]; mn=mp.mpf('9'); aty=None
M=6000
for i in range(1,M+1):
    y = y100*i/M; v = G(y,1)
    if v < mn: mn=v; aty=y
    if prev is not None and (prev[1]<0)!=(v<0): crossings.append((prev[0],y))
    prev=(y,v)
print(f"G(y,1) on (0,y100]: min={mp.nstr(mn,12)} at y={mp.nstr(aty,12)} (N={mp.nstr(mpp/(2*aty**2),8)}); sign crossings: {[(mp.nstr(mpp/(2*a**2),8), mp.nstr(mpp/(2*b**2),8)) for a,b in crossings]}")
# and G(y,1) at N=1e4, 1e6 for the sweep claim
for Nv in [400, 1e3, 1e4, 1e6]:
    y=mp.sqrt(mpp/(2*Nv)); print(f"   G(y(N={Nv}),1) = {mp.nstr(G(y,1),8)}")
