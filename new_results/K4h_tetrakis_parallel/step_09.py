step02 = r'''"""step_02.py -- K4h recon: outer splits, tanh-sinh weight cutoffs, node counts, cost calibration."""
import sys, os, json, time
sys.dont_write_bytecode = True
K4E = "<project>/research/08_theory/K4e_tetrakis_d3"
sys.path.insert(0, K4E)
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 60
from k4c_driver import outer_splits
from k4_quad import tanh_sinh_nodes

s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]

print("== outer splits for NA ==")
spl = outer_splits(NA)
splits = [mpf(0)] + spl + [mp.pi/4]
print("split points:", [mp.nstr(s,12) for s in splits])
lens = [splits[i+1]-splits[i] for i in range(len(splits)-1)]
print("n_intervals:", len(lens), " lengths:", [mp.nstr(L,12) for L in lens])
print("max outer interval length:", mp.nstr(max(lens),12))

# inner theta interval: [atan(1/sin phi), pi/2); thlo minimized at phi=pi/4
thlo_min = mp.atan(1/mp.sin(mp.pi/4))
inner_max = mp.pi/2 - thlo_min
print("thlo_min =", mp.nstr(thlo_min,15), " max inner interval length =", mp.nstr(inner_max,15))

def w_ts(t, a, b, h):
    """tanh-sinh weight at parameter t (same formula as k4_quad)."""
    sh = mp.sinh(t); ch = mp.cosh(t)
    cp = mp.cosh(mp.pi/2*sh)
    return (b-a)/2 * (mp.pi/2*ch)/(cp*cp) * h

def tstar(a, b, h, WMIN):
    """largest t with w(t) >= WMIN (bisection)."""
    lo, hi = mpf(0), mpf('20')
    for _ in range(300):
        mid = (lo+hi)/2
        if w_ts(mid, a, b, h) >= WMIN: lo = mid
        else: hi = mid
    return lo

print()
print("== weight cutoff analysis (t* such that w(t*) = WMIN) ==")
cases = [
    ("validation h=1/96  WMIN=1e-45 ", mpf(1)/96,  mpf('1e-45')),
    ("h384       h=1/384 WMIN=1e-125", mpf(1)/384, mpf('1e-125')),
    ("production h=1/768 WMIN=1e-215", mpf(1)/768, mpf('1e-215')),
]
for tag, h, W in cases:
    for (a,b,lab) in [(mpf(0), max(lens), "maxouter"), (mpf(0), inner_max, "maxinner"), (mpf(0), min(lens), "minouter")]:
        ts = tstar(a, b, h, W)
        print("  %s %s len=%s  t*=%s" % (tag, lab, mp.nstr(b-a,6), mp.nstr(ts,8)))
    print()

print("== recommended ranges T ==")
# production: need T >= t* for the largest interval (weights then reach WMIN before T cuts)
h768 = mpf(1)/768
ts_prod = max(tstar(mpf(0), max(lens), h768, mpf('1e-215')), tstar(mpf(0), inner_max, h768, mpf('1e-215')))
print("production t*_max =", mp.nstr(ts_prod,10), "-> choose T=6.0; check w(6.0) vs WMIN:")
for (a,b,lab) in [(mpf(0), max(lens), "maxouter"), (mpf(0), inner_max, "maxinner")]:
    print("   w(6.0) %s = %s (WMIN=1e-215)" % (lab, mp.nstr(w_ts(mpf('6.0'),a,b,h768),3)))
h96 = mpf(1)/96
ts_val = max(tstar(mpf(0), max(lens), h96, mpf('1e-45')), tstar(mpf(0), inner_max, h96, mpf('1e-45')))
print("validation t*_max =", mp.nstr(ts_val,10))
for T in ['4.4','4.5']:
    print("   w(%s) maxinner = %s (WMIN=1e-45)" % (T, mp.nstr(w_ts(mpf(T),mpf(0),inner_max,h96),3)))

print()
print("== outer node counts at production params (h=1/768, T=6.0, WMIN=1e-215) ==")
tot = 0
for i in range(len(lens)):
    nd = tanh_sinh_nodes(splits[i], splits[i+1], h768, mpf('6.0'), mpf('1e-215'))
    print("  interval %d len=%s nodes=%d" % (i, mp.nstr(lens[i],6), len(nd)))
    tot += len(nd)
print("  TOTAL outer nodes:", tot)

print()
print("== inner node count sample (typical phi) ==")
for phv in ['0.1','0.3','0.6']:
    ph = mpf(phv)
    thlo = mp.atan(1/mp.sin(ph))
    nd = tanh_sinh_nodes(thlo, mp.pi/2, h768, mpf('6.0'), mpf('1e-215'))
    print("  phi=%s thlo=%s innerlen=%s innernodes=%d" % (phv, mp.nstr(thlo,6), mp.nstr(mp.pi/2-thlo,6), len(nd)))

print()
print("== calibration from h384 run (22832 s, dps 160, T=5.4, WMIN=1e-125) ==")
tot384 = 0
for i in range(len(lens)):
    nd = tanh_sinh_nodes(splits[i], splits[i+1], mpf(1)/384, mpf('5.4'), mpf('1e-125'))
    tot384 += len(nd)
print("outer nodes h384:", tot384, " -> per-outer-node time = %.3f s" % (22832.0/tot384))

print()
print("== F_series evals/sec at dps 160 vs 270 (m=6, K=4) ==")
from k4_engine import FastModel
for dps in [160, 270]:
    mp.mp.dps = dps
    M = FastModel(NA, 4)
    ux,uy,uz = mpf('0.4'),mpf('0.5'),mp.sqrt(1-mpf('0.4')**2-mpf('0.5')**2)
    sig = (1,)*6
    M.F_series(ux,uy,uz,sig)
    n = 300
    t0=time.time()
    for _ in range(n): M.F_series(ux,uy,uz,sig)
    el=time.time()-t0
    print("  dps %d: %.1f evals/s (%.3f ms/eval)" % (dps, n/el, 1000*el/n))
mp.mp.dps = 60

# tanh-sinh node GENERATION cost (per outer node, inner list generation incl. dropped tails)
mp.mp.dps = 270
t0=time.time()
nd = tanh_sinh_nodes(mpf('1.0'), mp.pi/2, h768, mpf('6.0'), mpf('1e-215'))
el=time.time()-t0
print("  inner node-list generation at dps270: %.1f ms for %d kept nodes" % (1000*el, len(nd)))
print()
print("== rough single-process production estimate (to be refined by pilot) ==")
# per outer node: inner evals ~9k at dps270 + generation ~12k candidates
per_node = None
print("(see step_05 pilot)")
'''
open("step_02.py","w").write(step02)
import subprocess, sys
r = subprocess.run([sys.executable, "step_02.py"], capture_output=True, text=True)
print(r.stdout); print("ERR:", r.stderr[-3000:] if r.returncode else "(clean)")
