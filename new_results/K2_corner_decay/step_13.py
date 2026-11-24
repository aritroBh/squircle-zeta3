
import mpmath as mp
import time, json, os

mp.mp.dps = 140
DPS = mp.mp.dps
P   = 1024          # FFT points
KMAX = 300
R   = mp.mpf('0.9') # Cauchy radius
H   = mp.mpf('0.025')
JTS = 224           # tanh-sinh nodes per side  (t in [-5.6, 5.6])
M1, J1 = 6, 6

def gamma_coeffs(c, M1):
    Q = [mp.mpf(1), -2*c, mp.mpf(1)] + [mp.mpf(0)]*(M1-2)
    b = [mp.sqrt(Q[0])]
    for n in range(1, M1+1):
        s = Q[n] - mp.fsum(b[k]*b[n-k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b

# ---- global tanh-sinh grid for g in (0,1): g = (1+tanh(pi/2 sinh t))/2
t0 = time.time()
grid = []
for jj in range(-JTS, JTS+1):
    t = jj*H; st = mp.sinh(t)
    th = mp.tanh(mp.pi/2*st)
    g = (1+th)/2
    w = (mp.pi/4)*mp.cosh(t)/(mp.cosh(mp.pi/2*st)**2)
    if g <= 0 or g >= 1: continue
    grid.append((g, w, mp.log(g), mp.log1p(g)))
print("grid: %d nodes, built %.1fs" % (len(grid), time.time()-t0))

def compute_fvals(alpha):
    c = mp.cos(alpha)
    gam = gamma_coeffs(c, M1)
    f = []
    for p in range(P//2+1):
        u = R*mp.exp(2*mp.pi*mp.j*p/P)
        beta = [mp.mpc(1)]
        for j in range(J1):
            beta.append(beta[-1]*(-(1+u+j))/(j+1))
        ana = mp.mpc(0)
        for m in range(1, M1+1): ana += gam[m]/(m-(m-1)*u)
        for j in range(1, J1+1): ana += beta[j]/(u+j)
        tot = mp.mpc(0)
        for (g, w, lg, l1g) in grid:
            A = mp.exp(lg*(u-1))
            B = 1/A
            lC = -(1+u)*l1g
            C = mp.exp(lC)
            Q = 1 - 2*c*B + B*B
            sQ = mp.sqrt(Q)
            Cm1 = C*((Q-1)/(sQ+1)) + mp.expm1(lC)   # C*sqrt(Q)-1 stably
            S = Cm1*A
            Bp = mp.mpc(1)
            for m in range(1, M1+1):
                S -= gam[m]*Bp; Bp *= B
            Ap = A*g
            for j in range(1, J1+1):
                S -= beta[j]*Ap; Ap *= g
            tot += w*S
        f.append(tot*H + ana)
    return f

def dft(f):
    # a_k = (1/P)[ f0 + (-1)^k f_{P/2} + 2 Re sum_{p=1}^{P/2-1} f_p w^{-kp} ]
    half = P//2
    om = mp.exp(-2*mp.pi*mp.j/P)
    ak = []
    max_im = mp.mpf(0)
    for k in range(KMAX+1):
        s = mp.mpc(0)
        wk = om**k
        wp = mp.mpc(1)
        for p in range(1, half):
            wp *= wk
            s += f[p]*wp
        a = (f[0] + ((-1)**k)*f[half] + 2*s.real)/P
        max_im = max(max_im, abs(2*s.imag/P))
        ak.append(a * R**(-k))
    return ak, max_im

angles = [('pi_3', mp.pi/3), ('pi_2', mp.pi/2), ('2pi_3', 2*mp.pi/3),
          ('3pi_4', 3*mp.pi/4), ('5pi_6', 5*mp.pi/6)]

os.makedirs('out', exist_ok=True)
for name, alpha in angles:
    t0 = time.time()
    f = compute_fvals(alpha)
    ak, max_im = dft(f)
    dt = time.time()-t0
    with open('out/mk_%s.tsv' % name, 'w') as fh:
        fh.write("# M_k(alpha=%s*pi), Cauchy/FFT P=%d r=%s dps=%d\n" % (name, P, mp.nstr(R,3), DPS))
        for k, v in enumerate(ak):
            fh.write("%d\t%s\n" % (k, mp.nstr(v, 60)))
    print("%-7s done %.1fs  max|Im a_k|=%.1e   M0=%s" % (name, dt, max_im, mp.nstr(ak[0], 30)))
    print("   M50=%s\n   M100=%s\n   M200=%s\n   M300=%s" % tuple(mp.nstr(ak[k], 25) for k in (50,100,200,300)))
