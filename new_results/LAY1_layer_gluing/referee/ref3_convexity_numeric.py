# Referee item 3: independent numerical min of h + h'' for the soft-max, dps 300 (own code).
# h is computed from the raw vertex-sum definition; h'' by mpmath.diff (high-precision finite differences).
from mpmath import mp, mpf, cos, sin, exp, pi, diff, nstr, findroot
mp.dps = 300
V = [(1,1),(1,-1),(-1,1),(-1,-1)]
def h(th, tau):
    n = (cos(th), sin(th))
    d = [v[0]*n[0]+v[1]*n[1] for v in V]
    m = max(d)
    return m*sum(exp((x-m)/tau) for x in d)**tau
def k(th, tau):  # h + h''
    return h(th,tau) + diff(lambda t: h(t,tau), th, 2)
for tau in [mpf(1), mpf('0.5'), mpf('0.2'), mpf('0.1'), mpf('0.05'), mpf('0.02')]:
    N = 400
    grid = [pi/4*i/N for i in range(N+1)]
    vals = [k(t,tau) for t in grid]
    i = min(range(N+1), key=lambda j: vals[j])
    # refine with golden-section on neighbouring cell
    a = grid[max(i-1,0)]; b = grid[min(i+1,N)]
    gr = (mp.sqrt(5)-1)/2
    for _ in range(80):
        c = b - gr*(b-a); d = a + gr*(b-a)
        if k(c,tau) < k(d,tau): b = d
        else: a = c
    tm = (a+b)/2; km = k(tm,tau)
    # also check the extra points: theta=0 (kink) and theta=pi/4, and slightly negative theta (other side of kink)
    print(f"tau={nstr(tau,3):>5}  min(h+h'') on grid at theta={nstr(grid[i],6)}: {nstr(vals[i],8)} ; refined theta={nstr(tm,10)} -> {nstr(km,12)}"
          f" ; h+h''(0)={nstr(k(mpf(0),tau),8)}  h+h''(pi/4)={nstr(k(pi/4,tau),8)}  h+h''(-0.01)={nstr(k(mpf('-0.01'),tau),8)}")
