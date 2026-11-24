import sys, math, time
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, acb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')
iv, h3encl, pi, rt2 = st.iv, st.h3encl, st.pi, st.rt2
sq2 = rt2; a2 = pi*sq2/4
import mpmath as mp
mp.mp.dps = 60
r2 = mp.sqrt(2); mpp = mp.pi

cells = []; N0 = 100.0
while N0 < 400*(1-1e-12):
    N1 = min(N0*1.05, 400.0); cells.append((N0,N1)); N0 = N1

def coded_mS(y1,y2,s1,s2):
    y = iv(y1,y2); sig = iv(s1,s2)
    rho = 1 + y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)
def truth_pt(y, sig):
    y = arb(y); sig = arb(sig)
    rho = 1 + y*y*sig; t = rho*y; N = pi/(2*y*y)
    w = 1 - (1+sq2*t+t*t)*(-sq2*t).exp(); psi = -(1-w).log()
    return (N/2)*psi - (a2/t)*(rho-1)**2

def saddle_sig(yv):
    """sigma*(y): solve t/rho(t) = y, rho(t)=sqrt(1+t^2/(1+rt2 t)); return sigma, t, rho."""
    yv = mp.mpf(yv)
    f = lambda tt: tt/mp.sqrt(1+tt**2/(1+r2*tt)) - yv
    t = mp.findroot(f, yv*(1+yv**2/2))
    rho = mp.sqrt(1+t**2/(1+r2*t))
    return (rho-1)/yv**2, t, rho

# (a) saddle containment per cell
print("(a) saddle containment in subcell k=1 ([0.25,0.5]) mS box:")
sg_lo_all = mp.mpf(9); sg_hi_all = mp.mpf(0); w_lo = None; w_hi = None
for (n1,n2) in cells:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    mS = coded_mS(y1,y2,0.25,0.5); mlo = arb(mS.lower()); mup = arb(mS.upper())
    for yv in [y1, (y1+y2)/2, y2]:
        sg, tt, rho = saddle_sig(yv)
        sg_lo_all = min(sg_lo_all, sg); sg_hi_all = max(sg_hi_all, sg)
        assert 0.25 < sg < 0.5, f"saddle outside subcell k=1 at y={yv}"
        tv = truth_pt(yv, float(sg))
        m_lo = arb(tv.lower()) - mlo; m_hi = mup - arb(tv.upper())
        if w_lo is None or arb(m_lo.upper()) < arb(w_lo[0].upper()): w_lo = (m_lo,(n1,n2,yv,float(sg)))
        if w_hi is None or arb(m_hi.upper()) < arb(w_hi[0].upper()): w_hi = (m_hi,(n1,n2,yv,float(sg)))
print(f"   sigma* range over all 29 cells (y endpoints+mid): [{mp.nstr(sg_lo_all,12)}, {mp.nstr(sg_hi_all,12)}]  (R1 claim [0.428579581, 0.5))")
print(f"   WORST saddle lower margin = {float(arb(w_lo[0].mid())):+.6e}  at cell [{w_lo[1][0]:.3f},{w_lo[1][1]:.3f}], y={w_lo[1][2]:.7f}, sig*={w_lo[1][3]:.6f}")
print(f"   WORST saddle upper margin = {float(arb(w_hi[0].mid())):+.6e}  at cell [{w_hi[1][0]:.3f},{w_hi[1][1]:.3f}], y={w_hi[1][2]:.7f}, sig*={w_hi[1][3]:.6f}")

# detail for first cell (compare R5: margins +0.00522/+0.00884 for cell [100,105])
y1 = math.sqrt(math.pi/(2*105)); y2 = math.sqrt(math.pi/(2*100))
mS = coded_mS(y1,y2,0.25,0.5)
sg,tt,rho = saddle_sig(y2)
tv = truth_pt(y2, float(sg))
print(f"   cell [100,105] saddle: sig*={mp.nstr(sg,8)}, truth={float(arb(tv.mid())):.10f}, mS=[{float(arb(mS.lower())):.10f},{float(arb(mS.upper())):.10f}], margins +{float(arb(tv.lower())-arb(mS.lower())):.5f}/+{float(arb(mS.upper())-arb(tv.upper())):.5f}")
