
import mpmath as mp, json
mp.mp.dps = 60
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
mine = json.load(open('SRF1_n3_K60.json'))
tau = [mp.mpf(x) for x in h1['eta']]   # S/(4sqrt3) in tau, K=80
d   = [mp.mpf(x) for x in h1['d']]     # S/(4sqrt3) in d,   K=80
w   = [mp.mpf(x) for x in mine['w']]   # log[(2-p)S/(4sqrt3)] in tau, K=60 (mine)

A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2)   # predicted 0.62166...
print("predicted Borel |t_edge| =", mp.nstr(1/A3, 12), "  A = 2/sqrt(pi^2+log^2 2) =", mp.nstr(A3,12))

def analyze(a, name, Kmax):
    print(f"--- {name} ---")
    print(" k   |a_k|^(1/k)/k     (|a_k|/k!)^(1/k)   |a_{k+1}/a_k|/(k+1)")
    for k in range(10, Kmax, 10):
        rt = abs(a[k])**(1/k)/k
        bf = (abs(a[k])/mp.factorial(k))**(1/k)
        ra = abs(a[k+1]/a[k])/(k+1)
        print(f"{k:3d}  {mp.nstr(rt,10)}  {mp.nstr(bf,10)}  {mp.nstr(ra,10)}")
analyze(tau, "tau-series S/(4sqrt3), n=3", 79)
analyze(d,   "d-series S/(4sqrt3), n=3", 79)
analyze(w,   "log-series w_k, n=3", 59)
# Domb-Sykes slope: ratio |a_{k+1}/a_k| vs k should be linear with slope A
print()
print("Domb-Sykes: |a_{k+1}/a_k| = A*(k + 1 + beta) fit, last window k=70..79 (tau), k=50..59 (w)")
for a, nm, lo, hi in [(tau,'tau',70,79),(d,'d',70,79),(w,'w',50,59)]:
    ks = list(range(lo, hi))
    xs = [mp.mpf(k) for k in ks]; ys = [abs(a[k+1]/a[k]) for k in ks]
    npts = len(ks); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(x*y for x,y in zip(xs,ys))
    slope = (npts*sxy-sx*sy)/(npts*sxx-sx*sx); inter = (sy-slope*sx)/npts
    print(f"  {nm}: slope={mp.nstr(slope,8)} (pred {mp.nstr(A3,8)}), intercept={mp.nstr(inter,6)} -> beta={mp.nstr(inter/slope-1,6)}")
