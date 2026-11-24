"""RULQ script 3 (tasks B, C): optimal truncation of pi_p/4 ~ 1 + sum N_k p^{-k-1}.
Input: N_k (k <= 100, 110 digits) from ../PIPDIV_divergence_proof/pd1_N_K100_dps110_deg6_S760.txt.
 1. pi_p/4 by the split v-form (and the w-form as an independent check); pi_2 = pi.
 2. R_K(p) = pi_p/4 - 1 - sum_{k<K} N_k h^{k+1}; optimal K and smallest-term K.
 3. Models: lead   R^lead_K = 2Re[h^2 i^{1-K} int_0^oo e^{-pi p s} s^{K-2}/(1+i s) ds]
            simple R^simp_K = s_K h^2 Gamma(K-1) (pi p)^{1-K},  s_K = +1,-1,-1,+1 for K = 1,2,3,0 mod 4
            sing   R^sing_K = M(h) - sum_{k<K} S_k h^{k+1},  M(h) = 2Re[h^2 int_0^oo e^{-v/h}(i pi - v)^{-1-h} dv]
 4. The rigorous bound of Theorem B2 (result.md Sec. 5), minimised over its free parameters.
Usage: ru3_borel.py"""
import mpmath as mp
mp.mp.dps = 60
N = [mp.mpf(l) for l in open('../PIPDIV_divergence_proof/pd1_N_K100_dps110_deg6_S760.txt')]
KMAX = len(N) - 1

def L8_v(p, q):
    p = mp.mpf(p); q = mp.mpf(q); h = 1 / p; r = q * (1 - h)
    f = lambda v: mp.exp(-h * v) * (1 + mp.exp(-v)) ** (-1 - h) * mp.expm1(mp.log1p(mp.exp(-r * v)) / q)
    pts = sorted(set([0] + [mp.mpf(4) ** i / r for i in range(0, 6)] + [mp.mpf(1), 4, 16, 64, 256, 1024]))
    return mp.mpf(2) ** (-h) + h * mp.quad(f, pts)
def L8_w(p, q):
    p = mp.mpf(p); q = mp.mpf(q); e = q * (1 - 1 / p); W = mp.mpf(2) ** (-1 / p)
    f = lambda w: (1 + (w ** p / (1 - w ** p)) ** e) ** (1 / q)
    pts = [0, W / 2] + [W * (1 - mp.mpf(1) / 4 ** i) for i in range(1, 6)] + [W]
    return mp.quad(f, pts)

# singular coefficients S_k (lambda = 1): S_k = 2Re sum_m y0^{-1-m} [h^{k-2-m}] (e^{-h Log y0} (1+h)_m)
y0 = mp.mpc(0, mp.pi); w = -mp.log(y0)
n = KMAX + 1
eh = [w ** i / mp.factorial(i) for i in range(n)]
def smul(a, b, m):
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1)) for k in range(m)]
pochs = [[mp.mpf(1)] + [mp.mpf(0)] * (n - 1)]
for m in range(1, n): pochs.append(smul(pochs[-1], [mp.mpf(m), mp.mpf(1)] + [mp.mpf(0)] * (n - 2), n))
S = [mp.mpf(0)] * n
for k in range(2, n):
    s = mp.mpc(0)
    for m in range(0, k - 1):
        i = k - 2 - m
        s += y0 ** (-1 - m) * mp.fsum(eh[a] * pochs[m][i - a] for a in range(i + 1))
    S[k] = 2 * mp.re(s)
def Mfun(h):
    f = lambda v: mp.exp(-v / h) * mp.exp((-1 - h) * mp.log(y0 - v))
    return 2 * mp.re(h ** 2 * mp.quad(f, [0, h, 4 * h, 16 * h, 64 * h, mp.inf]))

def lead(p, K):
    p = mp.mpf(p); h = 1 / p
    f = lambda s: mp.exp(-mp.pi * p * s) * s ** (K - 2) / (1 + 1j * s)
    sstar = (K - 2) / (mp.pi * p)
    I = mp.quad(f, [0, sstar / 2, sstar, 2 * sstar, 4 * sstar, mp.inf])
    return 2 * mp.re(h ** 2 * mp.power(1j, 1 - K) * I)
def simple(p, K):
    sK = {1: 1, 2: -1, 3: -1, 0: 1}[K % 4]
    return sK * mp.gamma(K - 1) * (mp.pi * p) ** (1 - K) / mp.mpf(p) ** 2

# ---------------- rigorous bound (Theorem B2)
def cN(Nn):
    return mp.sqrt(mp.pi) * mp.gamma(mp.mpf(Nn - 1) / 2) / mp.gamma(mp.mpf(Nn) / 2)
def remfac(M, a, r, Nn, h):   # Lemma R: |Rem_N| <= M/pi (1+c_N) e^{a r} N! (h/(1-ah))^{N+1} r^{-N}
    return M / mp.pi * (1 + cN(Nn)) * mp.exp(a * r) * mp.factorial(Nn) * (h / (1 - a * h)) ** (Nn + 1) * r ** (-Nn)
def bound(p, K, rho, delta, rho2, delta2):
    h = 1 / mp.mpf(p)
    if not (h < rho and (1 + rho) * h < 1 and h < rho2 < 1): return mp.inf
    R0 = 2 * (h * mp.log(2)) ** (K + 1) / mp.factorial(K + 1)
    r1 = mp.pi - delta; s1 = mp.sqrt(2) * mp.sin(delta / 2); Lam = mp.pi + mp.log(2 / s1)
    M1 = mp.exp(rho * Lam) / s1; a1 = 1 + rho
    R1 = mp.fsum(h ** j * rho ** (-j) * remfac(M1, a1, r1, K - 2 - j, h) for j in range(0, K - 3))
    R1 += 2 * M1 * h * (h / rho) ** (K - 3) / ((1 - h / rho) * (1 - a1 * h))
    R1 *= h ** 2
    r2 = 2 * mp.pi - delta2; s2 = mp.sqrt(2) * mp.sin(delta2 / 4)
    M2 = rho2 ** 2 * mp.exp(rho2 * (1 + mp.pi)) * s2 ** (-1 - rho2) * (mp.zeta(2 - rho2) - 1); a2 = 1 + rho2 / 2
    R2 = mp.fsum(h ** j * rho2 ** (-j) * remfac(M2, a2, r2, K - j, h) for j in range(2, K - 1))
    R2 += 2 * M2 * h * (h / rho2) ** (K - 1) / ((1 - h / rho2) * (1 - a2 * h))
    return R0 + R1 + R2, R0, R1, R2
def best_bound(p, K):
    best = None
    h = 1 / mp.mpf(p)
    for rho in [h * t for t in (1.5, 2, 2.5, 3, 4, 5, 6)] + [0.6, 0.7, 0.8, 0.9]:
        for delta in [mp.mpf(c) / p for c in (0.25, 0.5, 1, 1.5, 2, 3)]:
            if delta >= mp.pi / 2: continue
            b = bound(p, K, mp.mpf(rho), delta, mp.mpf(0.9), mp.mpf(1))
            if b is mp.inf or b == mp.inf: continue
            if best is None or b[0] < best[0][0]: best = (b, rho, delta)
    # optimise the n >= 2 parameters at the best (rho, delta)
    b, rho, delta = best
    for rho2 in (0.6, 0.7, 0.8, 0.9, 0.95):
        for d2 in (0.5, 1, 1.5, 2):
            bb = bound(p, K, rho, delta, mp.mpf(rho2), mp.mpf(d2))
            if bb != mp.inf and bb[0] < best[0][0]: best = (bb, rho, delta)
    return best

if __name__ == '__main__':
    print('pi_p/4 values (v-form) and checks:')
    vals = {}
    for p in [2, 3, 4, 5, 6, 7, 8, 10, 12]:
        v = L8_v(p, p); vals[p] = v
        chk = mp.pi / 4 if p == 2 else L8_w(p, p)
        print('  p=%2d  pi_p/4 = %s   |v-w| (p=2: |v - pi/4|) = %s' % (p, mp.nstr(v, 45), mp.nstr(abs(v - chk), 3)))
    print('\nS_k check (S_k/N_k, k = 20, 40, 60):', [mp.nstr(S[k] / N[k], 10) for k in (20, 40, 60)])
    print('\nTable B (truncation of pi_p/4 - 1 ~ sum N_k p^{-k-1})')
    print(' p  K_s  |t_Ks|     R_{K_s}      K_opt  R_{K_opt}    sqrt2 p^-5/2 e^-pi p   lead(K_s)    simple(K_s)  sing(K_s)    bound(K_s)')
    rows = []
    for p in [2, 3, 5, 10]:
        h = 1 / mp.mpf(p)
        t = [N[k] * h ** (k + 1) for k in range(KMAX + 1)]
        R = [vals[p] - 1 - mp.fsum(t[:K]) for K in range(KMAX + 1)]
        Ks = min(range(2, KMAX), key=lambda k: abs(t[k]))
        Ko = min(range(2, KMAX), key=lambda k: abs(R[k]))
        ref = mp.sqrt(2) * mp.mpf(p) ** (-2.5) * mp.exp(-mp.pi * p)
        Mh = Mfun(h)
        sing = Mh - mp.fsum(S[k] * h ** (k + 1) for k in range(Ks))
        bb = best_bound(p, Ks)
        print('%2d  %3d  %.3e  % .4e   %3d  % .4e   %.4e        % .4e  % .4e  % .4e  %.3e' % (
            p, Ks, float(abs(t[Ks])), float(R[Ks]), Ko, float(R[Ko]), float(ref), float(lead(p, Ks)),
            float(simple(p, Ks)), float(sing), float(bb[0][0])))
        print('      bound parts R0=%.2e R1=%.2e R2=%.2e at rho=%s delta=%s' % (float(bb[0][1]), float(bb[0][2]), float(bb[0][3]),
              mp.nstr(bb[1], 3), mp.nstr(bb[2], 3)))
    print('\nTable C: R_K for K near pi p (p = 10, 12), ratio to models')
    for p in [8, 10, 12]:
        h = 1 / mp.mpf(p)
        t = [N[k] * h ** (k + 1) for k in range(KMAX + 1)]
        Mh = Mfun(h)
        for K in range(int(mp.pi * p) - 3, int(mp.pi * p) + 6):
            RK = vals[p] - 1 - mp.fsum(t[:K])
            sing = Mh - mp.fsum(S[k] * h ** (k + 1) for k in range(K))
            print('  p=%2d K=%2d  R_K=% .6e  R/lead=%.5f  R/simple=%.5f  R/sing=%.5f  |t_K|=%.3e' % (
                p, K, float(RK), float(RK / lead(p, K)), float(RK / simple(p, K)), float(RK / sing), float(abs(t[K]))))
    print('\nTable C2: R at K_s scaled by sqrt2 p^{-5/2} e^{-pi p} (p from 3 to 14)')
    for p2 in range(6, 29):
        p = mp.mpf(p2) / 2; h = 1 / p
        v = L8_v(p, p)
        t = [N[k] * h ** (k + 1) for k in range(KMAX + 1)]
        Ks = min(range(2, KMAX), key=lambda k: abs(t[k]))
        RK = v - 1 - mp.fsum(t[:Ks])
        ref = mp.sqrt(2) * p ** (-2.5) * mp.exp(-mp.pi * p)
        nu = Ks - 1 - mp.pi * p
        print('  p=%5s K_s=%2d nu=K_s-1-pi p=% .3f  R/ref=% .4f  s_K e^{nu^2/(2 pi p)}=% .4f  R/lead=%.4f' % (
            mp.nstr(p, 4), Ks, float(nu), float(RK / ref), float(simple(p, Ks) / ref), float(RK / lead(p, Ks))))
