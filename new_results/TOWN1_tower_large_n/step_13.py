import os, sys, time
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
sys.path.insert(0, WTS)
from fractions import Fraction as Fr
from math import comb, factorial
from fastred import QSI, ZERO, ONE, IQ, AQ, BQ

# ----- series arithmetic over QSI (series in s = v - i), truncated to order D (keep s^0..s^D) -----
def ser_mul(A, B, D):
    n = min(len(A)+len(B)-1, D+1)
    out = [ZERO]*n
    for a in range(len(A)):
        if a > D: break
        ca = A[a]
        if ca.is_zero(): continue
        lim = min(len(B), n-a)
        for b in range(lim):
            if not B[b].is_zero():
                out[a+b] = out[a+b] + ca*B[b]
    return out

def ser_pow_small(base, e, D):
    # base^e for e>=0, series, via repeated squaring
    res = [ONE] + [ZERO]*D
    if e == 0: return res
    b = list(base)
    while e:
        if e & 1: res = ser_mul(res, b, D)
        e >>= 1
        if e: b = ser_mul(b, b, D)
    return res

def ser_inv(A, D):
    a0 = A[0]
    a0i = a0.inv()
    out = [a0i] + [ZERO]*D
    for n in range(1, D+1):
        acc = ZERO
        for k in range(1, n+1):
            if k < len(A) and not A[k].is_zero() and not out[n-k].is_zero():
                acc = acc + A[k]*out[n-k]
        out[n] = (ZERO - acc)*a0i if not acc.is_zero() else ZERO
    return out

# Stirling numbers of 2nd kind, table up to n=maxn
def stirling_table(maxn):
    S = [[0]*(maxn+1) for _ in range(maxn+1)]
    S[0][0] = 1
    for n in range(1, maxn+1):
        for k in range(1, n+1):
            S[n][k] = k*S[n-1][k] + S[n-1][k-1]
    return S

def binom_half(j):
    r = Fr(1)
    for t in range(j):
        r *= (Fr(1,2) - t)
    return r / factorial(j)

# build the level-d series data
SQ2 = QSI(0,1,0,0)
def level_series(d, D, S2):
    # T = t(i+s), U = 1-T, W = q2(t(i+s)) = T^2+U^2
    # t(v) = (v^2+2v-1)/(4v); at v = i+s:
    # P(i+s) = (-2+2i) + (2+2i) s + s^2 ; 4(i+s) = 4i + 4s
    P = [QSI(-2,0,2,0), QSI(2,0,2,0), ONE] + [ZERO]*(D-2)   # coefficients of s^0, s^1, s^2
    P = P[:D+1]
    den = [QSI(0,0,4,0), QSI(4,0,0,0)] + [ZERO]*(D-1)
    den = den[:D+1]
    deni = ser_inv(den, D)
    T = ser_mul(P, deni, D)
    # U = 1 - T
    U = [ONE - T[0]] + [ZERO - c for c in T[1:]]
    # W = T^2 + U^2
    W = [a+b for a,b in zip(ser_mul(T,T,D), ser_mul(U,U,D))]
    # powers
    T2 = ser_mul(T, T, D)
    U2 = ser_mul(U, U, D)
    T2p = [[ONE]+[ZERO]*D]
    U2p = [[ONE]+[ZERO]*D]
    Wp  = [[ONE]+[ZERO]*D]
    for a in range(1, d+1):
        T2p.append(ser_mul(T2p[-1], T2, D))
        U2p.append(ser_mul(U2p[-1], U2, D))
        Wp.append(ser_mul(Wp[-1], W, D))
    # M[j][s] = T2p[j-s]*U2p[s]*Wp[d-j]  for j=1..d, s=0..j
    M = {}
    for j in range(1, d+1):
        Wdj = Wp[d-j]
        for s in range(0, j+1):
            M[(j,s)] = ser_mul(ser_mul(T2p[j-s], U2p[s], D), Wdj, D)
    # C_j = binom_half(j) * j! * 2^d
    C = [Fr(0)]*(d+1)
    for j in range(1, d+1):
        C[j] = binom_half(j) * factorial(j) * (2**d)
    # Common(s) = -(sqrt2/2) 8^d (i+s)^{2d-1} (2i+s)^{2-2d} / G(i+s), G = (v^2-a^2)(v^2-b^2)
    ips = [IQ, ONE] + [ZERO]*(D-1); ips = ips[:D+1]          # i + s
    v_pow = ser_pow_small(ips, 2*d-1, D)                      # (i+s)^{2d-1}
    # (2i+s)^{2-2d} = (2i)^{2-2d} (1 + s/(2i))^{2-2d}
    twoi = [QSI(0,0,2,0), ONE] + [ZERO]*(D-1); twoi = twoi[:D+1]
    twoi_inv = ser_inv(twoi, D)                               # (2i+s)^{-1}
    twoi_neg = ser_pow_small(twoi_inv, 2*d-2, D)              # (2i+s)^{2-2d}
    # G(i+s): G(v) = (v^2-a^2)(v^2-b^2) = v^4 - (a^2+b^2) v^2 + a^2 b^2
    # a^2 = (1+sqrt2)^2 = 3+2sqrt2 ; b^2 = 3-2sqrt2 ; a^2+b^2 = 6 ; a^2 b^2 = 1
    # G(v) = v^4 - 6 v^2 + 1
    # G(i+s): v^2 = (i+s)^2 = -1 + 2i s + s^2 ; v^4 = (v^2)^2
    v2 = ser_mul(ips, ips, D)
    v4 = ser_mul(v2, v2, D)
    G = [v4[t] - v2[t].scale(Fr(6)) + (ONE if t==0 else ZERO) for t in range(D+1)]
    Gi = ser_inv(G, D)
    Common = ser_mul(ser_mul(v_pow, twoi_neg, D), Gi, D)
    prefac = SQ2 * QSI(Fr(-(8**d), 2))
    Common = [c*prefac for c in Common]
    return dict(T=T, U=U, W=W, M=M, C=C, Common=Common, d=d, D=D)

def rho_key(LS, p, q, D, S2):
    d = LS['d']; M = LS['M']; C = LS['C']
    # Sp_ser = (1/(p! q!)) sum_j C_j sum_s S(p,s)S(q,j-s) M[(j,s)]
    pfact = factorial(p)*factorial(q)
    Sp = [ZERO]*(D+1)
    Spq = p+q  # = d
    for j in range(1, d+1):
        Cj = C[j]
        if Cj == 0: continue
        slo = max(0, j-q); shi = min(j, p)
        acc = [ZERO]*(D+1)
        for s in range(slo, shi+1):
            Sps = S2[p][s]; Sqj = S2[q][j-s]
            if Sps == 0 or Sqj == 0: continue
            coeff = Cj * Sps * Sqj
            Ms = M[(j,s)]
            for t in range(D+1):
                if not Ms[t].is_zero():
                    acc[t] = acc[t] + Ms[t].scale(coeff)
        for t in range(D+1):
            Sp[t] = Sp[t] + acc[t]
    inv_pf = Fr(1, pfact)
    Sp = [c.scale(inv_pf) for c in Sp]
    out = ser_mul(LS['Common'], Sp, D)
    if p > q:
        # fold: S_q(1-t): use M[(j,j-s)] with S(q,s)S(p,j-s)
        Sq = [ZERO]*(D+1)
        for j in range(1, d+1):
            Cj = C[j]
            if Cj == 0: continue
            slo = max(0, j-p); shi = min(j, q)
            acc = [ZERO]*(D+1)
            for s in range(slo, shi+1):
                Sqs = S2[q][s]; Spj = S2[p][j-s]
                if Sqs == 0 or Spj == 0: continue
                coeff = Cj * Sqs * Spj
                Ms = M[(j, j-s)]
                for t in range(D+1):
                    if not Ms[t].is_zero():
                        acc[t] = acc[t] + Ms[t].scale(coeff)
            for t in range(D+1):
                Sq[t] = Sq[t] + acc[t]
        Sq = [c.scale(inv_pf) for c in Sq]
        fold = ser_mul(LS['Common'], Sq, D)
        out = [a+b for a,b in zip(out, fold)]
    return out   # rho_{pq}(delta) = out[delta]

# validate for d=3 against towlib
S2 = stirling_table(20)
D = 4
t0=time.time()
LS = level_series(3, D, S2)
for (p,q) in [(3,0),(2,1)]:
    rho = rho_key(LS, p, q, D, S2)
    print(f"key ({p},{q}): rho(0..{D}) =")
    for delta in range(min(D, 2*(p+q)-4)+1):
        print("   delta", delta, rho[delta])
print("time", time.time()-t0)
