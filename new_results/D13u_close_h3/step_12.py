import sys, math, time
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, acb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')
iv, h3encl, pi, rt2 = st.iv, st.h3encl, st.pi, st.rt2
sq2 = rt2
t0 = time.time()

cells = []; N0 = 100.0
while N0 < 400*(1-1e-12):
    N1 = min(N0*1.05, 400.0); cells.append((N0,N1)); N0 = N1

def coded_mS(y1,y2,s1,s2):
    y = iv(y1,y2); sig = iv(s1,s2)
    rho = 1 + y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)

def truth_box(ybox, sbox):
    """interval outer-enclosure of RePhi(x0)+theta over the (y,sig) box, via exact identity,
       psi evaluated cancellation-free as -log(1-w), w = 1-(1+rt2 t+t^2)exp(-rt2 t)."""
    y = iv(*ybox) if isinstance(ybox, tuple) else arb(ybox)
    sig = iv(*sbox) if isinstance(sbox, tuple) else arb(sbox)
    rho = 1 + y*y*sig; t = rho*y; N = pi/(2*y*y)
    w = 1 - (1+sq2*t+t*t)*(-sq2*t).exp()
    psi = -(1-w).log()
    return (N/2)*psi - (pi*sq2/(4*t))*(rho-1)**2

NS, NY = 64, 64
glob_lo = None; glob_hi = None; nfail = 0; neval = 0
per_cell = []
for (n1,n2) in cells:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    cell_worst = None
    for k in range(4):
        s1, s2 = k/4, (k+1)/4
        mS = coded_mS(y1,y2,s1,s2); mlo = arb(mS.lower()); mup = arb(mS.upper())
        for i in range(NS):
            sa = s1+(s2-s1)*i/NS; sb = s1+(s2-s1)*(i+1)/NS
            for j in range(NY):
                ya = y1+(y2-y1)*j/NY; yb = y1+(y2-y1)*(j+1)/NY
                tv = truth_box((ya,yb),(sa,sb)); neval += 1
                m_lo = arb(tv.lower()) - mlo; m_hi = mup - arb(tv.upper())
                if arb(m_lo.upper()) < 0 or arb(m_hi.upper()) < 0: nfail += 1
                if glob_lo is None or arb(m_lo.upper()) < arb(glob_lo[0].upper()): glob_lo = (m_lo,(n1,n2,k,sa,sb,ya,yb))
                if glob_hi is None or arb(m_hi.upper()) < arb(glob_hi[0].upper()): glob_hi = (m_hi,(n1,n2,k,sa,sb,ya,yb))
                if cell_worst is None or arb(m_lo.upper()) < arb(cell_worst[0].upper()): cell_worst = (m_lo,k)
    per_cell.append((n1,n2,cell_worst))
print(f"interval-mesh certification: {neval} box evals ({NS}x{NS} per subcell),  ({time.time()-t0:.1f}s)")
print(f"boxes NOT contained: {nfail}")
print(f"WORST lower margin = [{arb(glob_lo[0].lower())}, {arb(glob_lo[0].upper())}]")
n1,n2,k,sa,sb,ya,yb = glob_lo[1]
print(f"   at cell [{n1:.4f},{n2:.4f}] subcell k={k}, sigma~[{sa:.5f},{sb:.5f}], y~[{ya:.7f},{yb:.7f}]")
print(f"WORST upper margin = [{arb(glob_hi[0].lower())}, {arb(glob_hi[0].upper())}]")
n1,n2,k,sa,sb,ya,yb = glob_hi[1]
print(f"   at cell [{n1:.4f},{n2:.4f}] subcell k={k}, sigma~[{sa:.5f},{sb:.5f}], y~[{ya:.7f},{yb:.7f}]")
print("\nper-cell worst lower margin (mid, and subcell k):")
for (n1,n2,(m,k)) in per_cell:
    print(f"  N [{n1:8.3f},{n2:8.3f}]  worst lower margin = {float(arb(m.mid())):+.6e}  (width {float(arb(m.upper())-arb(m.lower())):.2e}, subcell k={k})")
