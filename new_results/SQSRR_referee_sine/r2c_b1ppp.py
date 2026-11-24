# R2: b_n'(2), b_n''(2) straight from the DEFINITION of the squircle sine (paper, Sec. "Other Deformations"):
#   b_n(p) = 4/(n pi) int_0^1 cos(n Phi_p(y)) dy,  Phi_p(y) = (pi/2) F_p(y)/F_p(1) = (pi/2) I_{y^p}(1/p, 1-1/p)
# (regularized incomplete beta; F_p(y) = (1/p) B(y^p; 1/p, 1-1/p)).  Substitute y = sin(phi).
# p-derivatives: 7-point central finite differences on a FIXED node set (so the FD of the quadrature is the
# quadrature of the FD integrand), arb ball arithmetic at high precision.  Composite Gauss-Legendre, geometric
# grading at both ends.  Nothing from SQSR/AV8 is used.
import sys, time
from flint import arb, ctx
from mpmath import mp, mpf, legendre, cos as mcos, pi as mpi, nstr
P, G, HEXP, LEV = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
NS = [int(x) for x in sys.argv[5].split(',')]
ctx.prec = 600
mp.dps = 150
# Gauss-Legendre nodes on [-1,1]
xs, ws = [], []
for k in range(1, G+1):
    x = mcos(mpi*(k-mpf(1)/4)/(G+mpf(1)/2))
    for _ in range(100):
        p0, p1 = mpf(1), x
        for j in range(2, G+1): p0, p1 = p1, ((2*j-1)*x*p1-(j-1)*p0)/j
        dp = G*(x*p1-p0)/(x*x-1); dx = p1/dp; x -= dx
        if abs(dx) < mpf(10)**(-140): break
    p0, p1 = mpf(1), x
    for j in range(2, G+1): p0, p1 = p1, ((2*j-1)*x*p1-(j-1)*p0)/j
    dp = G*(x*p1-p0)/(x*x-1)
    xs.append(x); ws.append(2/((1-x*x)*dp*dp))
half = mpi/2; wdt = half/P
panels = [(i*wdt, (i+1)*wdt) for i in range(1, P-1)]
def graded(a, b, lev):   # grade toward a
    out = []; hi = b
    for _ in range(lev):
        lo = a + (hi-a)/2; out.append((lo, hi)); hi = lo
    out.append((a, hi)); return out
panels += graded(mpf(0), wdt, LEV)                                    # toward phi = 0
panels += [(half-b, half-a) for (a, b) in graded(mpf(0), wdt, LEV)]   # toward phi = pi/2 (stored as phi)
# nodes as u-distance from nearer endpoint to avoid cancellation
nodes = []
for (a, b) in panels:
    for x, w in zip(xs, ws):
        phi = (a+b)/2 + (b-a)/2*x
        nodes.append((phi, (b-a)/2*w))
print('nodes', len(nodes), flush=True)
h = arb(10)**(-HEXP)
offs = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
ps = [2 + k*h for k in offs]
PI = arb.pi(); HALF = PI/2
def phis(phi):
    """returns list over ps of (mode, value): mode 0 -> Phi itself, mode 1 -> pi/2 - Phi"""
    out = []
    if phi < mpi/4:
        y = arb(nstr(phi, 140)).sin()
        ly = y.log()
        for pp in ps:
            z = (pp*ly).exp()
            out.append((0, HALF*z.beta_lower(1/pp, 1-1/pp, 1)))
    else:
        u = arb(nstr(half-phi, 140))
        lc = (-2*(u/2).sin()**2).log1p()           # log cos u
        for pp in ps:
            z = -(pp*lc).expm1()                     # 1 - y^p
            out.append((1, HALF*z.beta_lower(1-1/pp, 1/pp, 1)))
    return out
t0 = time.time()
acc = {n: [arb(0)]*9 for n in NS}
for idx, (phi, w) in enumerate(nodes):
    ph = phis(phi)
    cphi = arb(nstr(mcos(phi), 140)); W = arb(nstr(w, 140))*cphi
    for n in NS:
        sg = 1 if (n-1)//2 % 2 == 0 else -1
        row = acc[n]
        for j, (mode, v) in enumerate(ph):
            c = (n*v).cos() if mode == 0 else sg*(n*v).sin()
            row[j] = row[j] + W*c
    if idx % 4000 == 0: print(idx, '%.0fs' % (time.time()-t0), flush=True)
# 9-point central stencils, O(h^8) for d1,d2 and O(h^6) for d3
c1 = [3, -32, 168, -672, 0, 672, -168, 32, -3]; c2 = [-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9]
c3 = [-7, 72, -338, 488, 0, -488, 338, -72, 7]
for n in NS:
    b = [4*v/(n*PI) for v in acc[n]]
    d1 = sum(c*v for c, v in zip(c1, b))/(840*h)
    d2 = sum(c*v for c, v in zip(c2, b))/(5040*h*h)
    d3 = sum(c*v for c, v in zip(c3, b))/(240*h**3)
    print('n', n, 'b(2)', b[4].str(40, radius=False))
    print('   b1 ', d1.str(60, radius=False))
    print('   b2 ', d2.str(60, radius=False))
    print('   b3 ', d3.str(60, radius=False), flush=True)
print('time %.0fs' % (time.time()-t0))
