import sys, math, time, subprocess
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')
iv, h3encl, pi, rt2 = st.iv, st.h3encl, st.pi, st.rt2
sq2 = rt2; a2 = pi*sq2/4

def coded_mS(y1,y2,s1,s2):
    y = iv(y1,y2); sig = iv(s1,s2)
    rho = 1 + y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)
def psi_f(t):
    w = 1 - (1+sq2*t+t*t)*(-sq2*t).exp(); return -(1-w).log()
def psip_f(t): return sq2*t*t/(1+sq2*t+t*t)
def truth_pt(y, sig):
    y = arb(y); sig = arb(sig)
    rho = 1 + y*y*sig; t = rho*y; N = pi/(2*y*y)
    return (N/2)*psi_f(t) - (a2/t)*(rho-1)**2
def dtruth_dy(y, sig):
    y = arb(y); sig = arb(sig); rho = 1+y*y*sig; t = rho*y; N = pi/(2*y*y)
    return (-(N/y)*psi_f(t) + (N/2)*psip_f(t)*(3*rho-2)
            - a2*(rho-1)*(4*y*sig*t - (rho-1)*(3*rho-2))/t**2)
def dtruth_ds(y, sig):
    y = arb(y); sig = arb(sig); rho = 1+y*y*sig; t = rho*y; N = pi/(2*y*y)
    return y**3*((N/2)*psip_f(t) - a2*(rho*rho-1)/t**2)
def truth_box_mv(ybox, sbox):
    yb = iv(*ybox); sb = iv(*sbox); ym = arb(yb.mid()); sm = arb(sb.mid())
    return truth_pt(ym, sm) + dtruth_dy(yb,sb)*(yb-ym) + dtruth_ds(yb,sb)*(sb-sm)

# (b) refine worst cell [392.0129138, 400], subcell k=3, at 256x256 + enclosure width diagnostic
n1,n2 = 392.0129138458655, 400.0
y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
mS = coded_mS(y1,y2,0.75,1.0); mlo = arb(mS.lower()); mup = arb(mS.upper())
NS=NY=256; w_lo=None; w_hi=None; maxw = arb(0)
for i in range(NS):
    sa = 0.75+0.25*i/NS; sb = 0.75+0.25*(i+1)/NS
    for j in range(NY):
        ya = y1+(y2-y1)*j/NY; yb = y1+(y2-y1)*(j+1)/NY
        tv = truth_box_mv((ya,yb),(sa,sb))
        wdt = arb(tv.upper())-arb(tv.lower())
        if arb(wdt.upper()) > arb(maxw.upper()): maxw = wdt
        m_lo = arb(tv.lower())-mlo; m_hi = mup-arb(tv.upper())
        if w_lo is None or arb(m_lo.upper()) < arb(w_lo.upper()): w_lo = m_lo
        if w_hi is None or arb(m_hi.upper()) < arb(w_hi.upper()): w_hi = m_hi
print(f"(b) cell [{n1},{n2}] k=3, 256x256 mean-value mesh:")
print(f"    certified worst lower margin = [{arb(w_lo.lower())}, {arb(w_lo.upper())}]")
print(f"    certified worst upper margin = [{arb(w_hi.lower())}, {arb(w_hi.upper())}]")
print(f"    max truth-enclosure width over mesh = {arb(maxw.upper())}")
# point value at the exact corner (y1, sig=1) for comparison
tv = truth_pt(y1, 1.0)
print(f"    point margin at (y1={y1:.8f}, sig=1): truth-mS.lower = {arb(arb(tv.lower())-mlo)}")
