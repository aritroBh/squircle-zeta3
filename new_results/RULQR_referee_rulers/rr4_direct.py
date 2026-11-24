"""RULQR script 4: (a) L_q(p) by direct quadrature of the geometric definition, x = cos^{2/p} t, y = sin^{2/p} t,
L_q = 8 int_0^{pi/4} (|x'|^q + |y'|^q)^{1/q} dt (symmetry t -> pi/2 - t; tanh-sinh handles the t^{2/p-1} endpoint),
against the Arb w-form of rr1 and, for q = p' = p/(p-1), against 8 Gamma(1+1/p)^2/Gamma(1+2/p);
(b) Table A'-type check: |L/8 - optimally truncated series| with the N_k of rr2;
(c) pi_p at p = 5, 10, 20 (Arb, 60 digits) and the optimally truncated sums (N_k: rr2 for k <= 50; PIPDIVR ball-arithmetic
N_k for p = 20, k up to 150)."""
import mpmath as mp
from flint import acb, arb, ctx
mp.mp.dps = 40; ctx.dps = 80

def L_direct(p, q):
    """x = cos^{2/p} t, y = sin^{2/p} t on [0, pi/4], times 8; t = u^{p/2} removes the t^{2/p-1} endpoint singularity."""
    p, q = mp.mpf(p), mp.mpf(q)
    def sp(u):
        if u == 0: return mp.mpf(1)   # limit of the integrand at u = 0
        t = u ** (p / 2); c, s = mp.cos(t), mp.sin(t)
        xd = 2 / p * c ** (2 / p - 1) * s; yd = 2 / p * s ** (2 / p - 1) * c
        return (xd ** q + yd ** q) ** (1 / q) * (p / 2) * u ** (p / 2 - 1)
    U = (mp.pi / 4) ** (2 / p)
    pts = [0, U / 2] + [U * (1 - mp.mpf(2) ** (-e)) for e in range(2, 40)] + [U]   # crossover at t = pi/4 is sharp for large q
    return 8 * mp.quad(sp, pts, maxdegree=10)

def L_graph(p, q):
    """second direct form: graph x = (1-y^p)^{1/p}, 0 <= y <= 2^{-1/p}, times 8."""
    p, q = mp.mpf(p), mp.mpf(q)
    f = lambda y: (1 + (y ** (p - 1) * (1 - y ** p) ** (1 / p - 1)) ** q) ** (1 / q)
    Y = mp.mpf(2) ** (-1 / p)
    return 8 * mp.quad(f, [0, Y / 2] + [Y * (1 - mp.mpf(2) ** (-e)) for e in range(2, 40)] + [Y], maxdegree=10)

def L_wform(p, q, dps=60):
    h = 1 / arb(p); q = arb(q); r = q * (1 - h); eta = 1 / q
    def f(w, a):
        A = (1 + (-w / r).exp()).log(analytic=a); B = (1 + (-w).exp()).log(analytic=a)
        return (-h * w / r - (1 + h) * A).exp() * (eta * B).expm1()
    W = int((dps + 20) * 2.31) + 10
    I = acb.integral(f, 0, W).real
    return 8 * ((arb(2) ** (-h)) + h / r * I)   # tail <= e^{-W} relative, negligible at this precision

print('(a) direct quadrature (40 dps; cos^{2/p} parametrisation and graph form) vs w-form (Arb) vs closed form for the dual ruler')
cases = [('3', '2'), ('5', '5'), ('4', '8'), ('10', '5'), ('6', '7'), ('2.5', '2.4'), ('4', '4/3'), ('1.5', '3'), ('7', '7/6'), ('20', '40')]
for ps, qs in cases:
    p = mp.mpf(ps); q = mp.mpf(qs) if '/' not in qs else None
    if '/' in qs: q = mp.mpf(qs.split('/')[0]) / mp.mpf(qs.split('/')[1])
    d = L_direct(p, q); g = L_graph(p, q)
    w = L_wform(ps, arb(qs.split('/')[0]) / arb(qs.split('/')[1]) if '/' in qs else qs)
    wm = mp.mpf(w.mid().str(45, radius=False))
    line = 'p=%-4s q=%-5s direct %s  |direct - wform| %.1e  |graph - wform| %.1e' % (ps, qs, mp.nstr(d, 30), float(abs(d - wm)), float(abs(g - wm)))
    if abs(q - p / (p - 1)) < 1e-30:
        cf = 8 * mp.gamma(1 + 1 / p) ** 2 / mp.gamma(1 + 2 / p)
        line += '  |direct - 8G(1+1/p)^2/G(1+2/p)| %.1e  area*2 = %s' % (float(abs(d - cf)), mp.nstr(cf, 12))
    print(line, flush=True)

rd = lambda f: [mp.mpf(l.split()[0]) for l in open(f) if l.strip()]
print('\n(b) |L/8 - S_{K*}| (S_{K*} = sum up to the smallest term), N_k from rr2 (k <= 50)')
for fam, qf, lam, plist in [('lam2', lambda p: 2 * p, 2, [4, 6]), ('lamhalf', lambda p: p / 2, mp.mpf(1) / 2, [8, 16]),
                            ('p1', lambda p: p + 1, 1, [8, 12]), ('gen', lambda p: p + mp.mpf(1) / 3 + 1 / p, 1, [8, 12])]:
    N = rd('rr2_N_%s_a.txt' % fam)
    for p in plist:
        p = mp.mpf(p); h = 1 / p
        q = qf(p)
        Lw = L_wform(str(p), arb(mp.nstr(q, 50)), 60)
        F = mp.mpf(Lw.mid().str(50, radius=False)) / 8
        t = [N[k] * h ** (k + 1) for k in range(len(N))]
        Ks = min(range(2, len(N)), key=lambda k: abs(t[k]))
        S = 1 + mp.fsum(t[:Ks])
        print('  %-8s p=%-3s K*=%2d  smallest term %.2e  |L/8 - S| %.2e  e^{-pi lam p} %.2e' % (fam, mp.nstr(p, 3), Ks, float(abs(t[Ks])),
              float(abs(F - S)), float(mp.exp(-mp.pi * lam * p))), flush=True)

print('\n(c) pi_p/4 (Arb w-form), optimal truncation; ref = sqrt2 p^{-5/2} e^{-pi p}')
Nr = rd('rr2_N_p_a.txt')
NP = [mp.mpf(l.split('[')[1].split('+/-')[0]) for l in open('../PIPDIVR_referee_divergence/rf1_N_K151_p700_S400.txt')]
mp.mp.dps = 60
for p in [5, 10, 20]:
    N = Nr if p < 20 else NP
    P4 = L_wform(str(p), str(p), 70) / 8        # L_p = 2 pi_p, so pi_p/4 = L_p/8
    P4 = mp.mpf(P4.mid().str(60, radius=False))
    h = mp.mpf(1) / p
    t = [N[k] * h ** (k + 1) for k in range(len(N))]
    Ks = min(range(2, len(N) - 1), key=lambda k: abs(t[k]))
    RK = lambda K: P4 - 1 - mp.fsum(t[:K])
    Kopt = min(range(2, min(len(N), 3 * Ks)), key=lambda K: abs(RK(K)))
    ref = mp.sqrt(2) * mp.mpf(p) ** (-2.5) * mp.exp(-mp.pi * p)
    print(' p=%2d pi_p/4 = %s' % (p, mp.nstr(P4, 40)))
    print('      K_s=%d |t_Ks|=%.3e R_Ks=%.4e  K_opt=%d R_Kopt=%.4e  ref=%.3e  |R_Ks|/ref=%.3f' % (Ks, float(abs(t[Ks])), float(RK(Ks)),
          Kopt, float(RK(Kopt)), float(ref), float(abs(RK(Ks)) / ref)))
    row = []
    for K in range(Ks - 4, Ks + 5):
        pred = mp.sqrt(2) * mp.cos(mp.pi * (K - 1) / 2 + mp.pi / 4) * mp.gamma(K - 1) / (mp.mpf(p) ** 2 * (mp.pi * p) ** (K - 1))
        row.append('K=%d R=%.3e R/Rsimple=%.3f' % (K, float(RK(K)), float(RK(K) / pred)))
    print('      ' + '\n      '.join(row), flush=True)
