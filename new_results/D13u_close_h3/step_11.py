import sys, math, random, time
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
from flint import arb, acb, ctx
ctx.dps = 70
import importlib
st = importlib.import_module('step_38')     # fixed certificate code (defs only; __main__ guarded)
iv, h3encl, pi, rt2 = st.iv, st.h3encl, st.pi, st.rt2
sq2 = rt2
t0 = time.time()

# ---- rebuild step_18's ratio-1.05 cover of [100,400] ----
cells = []; N0 = 100.0
while N0 < 400*(1-1e-12):
    N1 = min(N0*1.05, 400.0); cells.append((N0,N1)); N0 = N1
print(f"{len(cells)} cells (task says 29)")

def coded_mS(y1,y2,s1,s2):
    """EXACTLY step_38.py lines 50-59: box for RePhi(x0)+theta over (y,sig) box."""
    y = iv(y1,y2); sig = iv(s1,s2)
    rho = 1 + y*y*sig; y0 = rho*y; N = pi/(2*y*y)
    return (N*y0**3/2)*h3encl(y0)

def truth_identity(y, sig):
    """RePhi(x0)+theta = -(N/2)S(t) - (pi rt2/(4t))(rho-1)^2 = (N/2)psi(t) - ... ; arb (point or box)."""
    y = arb(y); sig = arb(sig)
    rho = 1 + y*y*sig; t = rho*y; N = pi/(2*y*y)
    # psi(t) = rt2 t - log(1+rt2 t + t^2), cancellation-free: psi = -log(1-w), w = 1-(1+rt2 t+t^2)e^{-rt2 t}
    w = 1 - (1+sq2*t+t*t)*(-sq2*t).exp()
    psi = -(1-w).log()
    return (N/2)*psi - (pi*sq2/(4*t))*(rho-1)**2

def truth_direct(y, sig):
    """direct from Phi (acb): Re[(-1+i) a x0 + N log(b x0/(1+b x0))] + theta."""
    y = arb(y); sig = arb(sig)
    rho = 1 + y*y*sig; x0 = 1/(rho*y); N = pi/(2*y*y); th = pi/(sq2*y)
    a = pi/(2*sq2); b = (acb(0,1)*pi/4).exp(); I1 = acb(0,1)
    x0c = acb(x0)
    return arb(((-1+I1)*a*x0c + N*(b*x0c/(1+b*x0c)).log()).real) + th

# ---- identity residual spot check (point evals) ----
random.seed(20240513)
worst_id = arb(0)
for _ in range(600):
    Nv = 100*4**random.random(); yv = math.sqrt(math.pi/(2*Nv)); sg = random.random()
    d = truth_identity(yv, sg) - truth_direct(yv, sg)
    r = abs(d)/abs(truth_direct(yv,sg))
    if arb(r.upper()) > arb(worst_id.upper()): worst_id = r
print(f"identity residual (600 random pts, arb/acb @70dps): worst rel = {arb(worst_id.upper())}")

# ---- main point check: all cells x 4 subcells; 5 fixed + 200 random sigma; y endpoints + 3 random y ----
worst_lo = None; worst_hi = None   # (margin, where)
npts = 0
for (n1,n2) in cells:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    for k in range(4):
        s1, s2 = k/4, (k+1)/4
        mS = coded_mS(y1,y2,s1,s2)
        mlo = arb(mS.lower()); mup = arb(mS.upper())
        sigs = [s1 + (s2-s1)*j/4 for j in range(5)] + [s1+(s2-s1)*random.random() for _ in range(200)]
        for sg in sigs:
            ys_ = [y1, y2] + [y1+(y2-y1)*random.random() for _ in range(3)]
            for yv in ys_:
                tv = truth_identity(yv, sg)
                m_lo = arb(tv.lower()) - mlo       # margin vs lower endpoint
                m_hi = mup - arb(tv.upper())       # margin vs upper endpoint
                npts += 1
                if worst_lo is None or arb(m_lo.upper()) < arb(worst_lo[0].upper()): worst_lo = (m_lo,(n1,n2,k,sg,yv))
                if worst_hi is None or arb(m_hi.upper()) < arb(worst_hi[0].upper()): worst_hi = (m_hi,(n1,n2,k,sg,yv))
print(f"point checks: {npts} evals  ({time.time()-t0:.1f}s)")
print(f"WORST lower margin (truth - mS.lower) = {arb(worst_lo[0].lower())} .. {arb(worst_lo[0].upper())}")
print(f"   at cell N in [{worst_lo[1][0]},{worst_lo[1][1]}], subcell k={worst_lo[1][2]} [{worst_lo[1][2]/4},{(worst_lo[1][2]+1)/4}], sigma={worst_lo[1][3]:.6f}, y={worst_lo[1][4]:.8f}")
print(f"WORST upper margin (mS.upper - truth) = {arb(worst_hi[0].lower())} .. {arb(worst_hi[0].upper())}")
print(f"   at cell N in [{worst_hi[1][0]},{worst_hi[1][1]}], subcell k={worst_hi[1][2]} [{worst_hi[1][2]/4},{(worst_hi[1][2]+1)/4}], sigma={worst_hi[1][3]:.6f}, y={worst_hi[1][4]:.8f}")
