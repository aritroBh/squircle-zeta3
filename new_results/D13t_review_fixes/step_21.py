import math, json, time
import mpmath as mp
from flint import arb, acb, ctx
mp.mp.dps = 60
ctx.dps = 70
pi_a = arb.pi(); rt2_a = arb(2).sqrt()
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)

# verbatim logic of fixed/step_38.py lines 168-177 (A1m box over the cell)
def A1m_box(y1, y2, s1, s2):
    y = iv(y1, y2); sig = iv(s1, s2)
    rho = 1 + y*y*sig
    zz = rt2_a*rho*y/(2+rt2_a*rho*y)
    zz_lo = arb(zz.lower()); zz_up = arb(zz.upper())
    fA1_lo = arb((1 - zz_lo.atan()/zz_lo).lower())
    fA1_up = arb((1 - zz_up.atan()/zz_up).upper())
    A1m = iv(fA1_lo/arb(y.upper()), fA1_up/arb(y.lower()))
    return (float(arb(A1m.lower())), float(arb(A1m.upper())),
            float(arb(zz_lo.lower())), float(arb(zz_up.upper())))

mpr2 = mp.sqrt(2)
def A1_true(yv, sgv):
    rho = 1 + yv*yv*sgv
    zz = mpr2*rho*yv/(2+mpr2*rho*yv)
    return (1 - mp.atan(zz)/zz)/yv, zz

# cells: near N=100, 1891 (sweep-worst), 1e4, 1e6, 1e9, 1e10-eps; mix of wide & sweep-narrow
cellsN = [(100.,105.),(100.,150.),(1890.7128907662395,1891.0910333443928),(1889.,1892.),
          (1e4,1.02e4),(1e6,1.0002e6),(1e8,1.0002e8),(4e9,1e10-1),(100.,400.)]
t0=time.time()
worst_viol = mp.mpf('-9'); worst_at=None; npts=0
worst_zz = mp.mpf('-9')
for (n1,n2) in cellsN:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    NS = 4
    for kk in range(NS):
        s1, s2 = kk/NS, (kk+1)/NS
        blo, bhi, zlo, zhi = A1m_box(y1*(1-1e-13), y2*(1+1e-13), s1, s2)
        G = 33
        for iy in range(G):
            yv = mp.mpf(y1) + (mp.mpf(y2)-mp.mpf(y1))*iy/(G-1)
            for isg in range(G):
                sgv = mp.mpf(s1) + (mp.mpf(s2)-mp.mpf(s1))*isg/(G-1)
                a1, zz = A1_true(yv, sgv)
                npts += 1
                viol = max(mp.mpf(blo)-a1, a1-mp.mpf(bhi))
                if viol > worst_viol: worst_viol = viol; worst_at=(n1,n2,s1,s2,float(yv),float(sgv))
                vzz = max(mp.mpf(zlo)-zz, zz-mp.mpf(zhi))
                worst_zz = max(worst_zz, vzz)
print(f"R1: {npts} sample points over {len(cellsN)} cells x 4 sigma-subcells")
print(f"  worst A1 containment violation = {mp.nstr(worst_viol,8)}  (<=0 => contained) at {worst_at}")
print(f"  worst zz-box violation         = {mp.nstr(worst_zz,8)}")
print(f"  ({time.time()-t0:.1f}s)")
