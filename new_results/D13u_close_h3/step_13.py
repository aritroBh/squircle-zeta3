import sys, math, time
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, acb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')
iv, h3encl, pi, rt2 = st.iv, st.h3encl, st.pi, st.rt2
sq2 = rt2; a2 = pi*sq2/4
t0 = time.time()

cells = []; N0 = 100.0
while N0 < 400*(1-1e-12):
    N1 = min(N0*1.05, 400.0); cells.append((N0,N1)); N0 = N1

def coded_mS(y1,y2,s1,s2):
    y = iv(y1,y2); sig = iv(s1,s2)
    rho = 1 + y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)

def psi_f(t):   # cancellation-free psi
    w = 1 - (1+sq2*t+t*t)*(-sq2*t).exp()
    return -(1-w).log()
def psip_f(t):  # psi'(t) = rt2 t^2/(1+rt2 t+t^2), no cancellation
    return sq2*t*t/(1+sq2*t+t*t)
def truth_pt(y, sig):
    y = arb(y); sig = arb(sig)
    rho = 1 + y*y*sig; t = rho*y; N = pi/(2*y*y)
    return (N/2)*psi_f(t) - (a2/t)*(rho-1)**2
def dtruth_dy(y, sig):
    y = arb(y); sig = arb(sig)
    rho = 1+y*y*sig; t = rho*y; N = pi/(2*y*y)
    return (-(N/y)*psi_f(t) + (N/2)*psip_f(t)*(3*rho-2)
            - a2*(rho-1)*(4*y*sig*t - (rho-1)*(3*rho-2))/t**2)
def dtruth_ds(y, sig):
    y = arb(y); sig = arb(sig)
    rho = 1+y*y*sig; t = rho*y; N = pi/(2*y*y)
    return y**3*((N/2)*psip_f(t) - a2*(rho*rho-1)/t**2)

# --- validate derivative formulas vs central finite differences (point) ---
import random; random.seed(7)
worst = arb(0)
for _ in range(300):
    Nv = 100*4**random.random(); yv = arb(math.sqrt(math.pi/(2*Nv))); sg = arb(random.random())
    h = arb(10)**(-20)
    fy = (truth_pt(yv+h,sg)-truth_pt(yv-h,sg))/(2*h); fs = (truth_pt(yv,sg+h)-truth_pt(yv,sg-h))/(2*h)
    r = max(abs(dtruth_dy(yv,sg)-fy)/max(abs(fy),arb(10)**(-30)), abs(dtruth_ds(yv,sg)-fs)/max(abs(fs),arb(10)**(-30)))
    if arb(r.upper()) > arb(worst.upper()): worst = r
print(f"derivative formulas vs finite differences: worst rel resid = {arb(worst.upper())}  (h=1e-20, so ~1e-14 floor expected)")

def truth_box_mv(ybox, sbox):
    """mean-value (centered) enclosure of truth over the (y,sig) box - rigorous."""
    yb = iv(*ybox); sb = iv(*sbox)
    ym = arb(yb.mid()); sm = arb(sb.mid())
    Dy = dtruth_dy(yb, sb); Ds = dtruth_ds(yb, sb)
    return truth_pt(ym, sm) + Dy*(yb-ym) + Ds*(sb-sm)

NS, NY = 64, 64
glob_lo = None; glob_hi = None; nfail = 0; neval = 0; per_cell=[]
for (n1,n2) in cells:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    cw = None
    for k in range(4):
        s1, s2 = k/4, (k+1)/4
        mS = coded_mS(y1,y2,s1,s2); mlo = arb(mS.lower()); mup = arb(mS.upper())
        for i in range(NS):
            sa = s1+(s2-s1)*i/NS; sb = s1+(s2-s1)*(i+1)/NS
            for j in range(NY):
                ya = y1+(y2-y1)*j/NY; yb = y1+(y2-y1)*(j+1)/NY
                tv = truth_box_mv((ya,yb),(sa,sb)); neval += 1
                m_lo = arb(tv.lower()) - mlo; m_hi = mup - arb(tv.upper())
                if arb(m_lo.upper()) < 0 or arb(m_hi.upper()) < 0: nfail += 1
                if glob_lo is None or arb(m_lo.upper()) < arb(glob_lo[0].upper()): glob_lo=(m_lo,(n1,n2,k,sa,sb,ya,yb))
                if glob_hi is None or arb(m_hi.upper()) < arb(glob_hi[0].upper()): glob_hi=(m_hi,(n1,n2,k,sa,sb,ya,yb))
                if cw is None or arb(m_lo.upper()) < arb(cw[0].upper()): cw=(m_lo,k)
    per_cell.append((n1,n2,cw))
print(f"mean-value mesh: {neval} box evals ({time.time()-t0:.1f}s); boxes not contained: {nfail}")
print(f"WORST lower margin = [{arb(glob_lo[0].lower())}, {arb(glob_lo[0].upper())}]")
n1,n2,k,sa,sb,ya,yb = glob_lo[1]; print(f"   cell [{n1:.4f},{n2:.4f}] k={k} sig~[{sa:.5f},{sb:.5f}] y~[{ya:.7f},{yb:.7f}]")
print(f"WORST upper margin = [{arb(glob_hi[0].lower())}, {arb(glob_hi[0].upper())}]")
n1,n2,k,sa,sb,ya,yb = glob_hi[1]; print(f"   cell [{n1:.4f},{n2:.4f}] k={k} sig~[{sa:.5f},{sb:.5f}] y~[{ya:.7f},{yb:.7f}]")
print("\nper-cell worst lower margin:")
for (n1,n2,(m,k)) in per_cell:
    print(f"  N [{n1:8.3f},{n2:8.3f}]  {float(arb(m.mid())):+.6e}  (encl width {float(arb(m.upper())-arb(m.lower())):.1e}, k={k})")
