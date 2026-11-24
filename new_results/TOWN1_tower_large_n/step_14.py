import os, sys, time
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
# progress log
with open('progress.log','a') as f:
    f.write(f"step_01 prototype rho series engine: d=3 matches recorded tower (delta 0..2)\n")

# Write the trunc_engine.py module
engine = r'''
"""trunc_engine: truncated tower engine at v = i, depth <= D, exact in Q(sqrt2, i).

Computes rho_m (initial data, top D+1 Laurent coefficients at v=i) directly from
Laurent/Taylor expansions (series in s = v - i), NOT from the full integrand PF.
Then applies the T1 depth-form recursion (towlib.Lambda_ij) to build the tower.

Theory: TOW1 proof_tower.md sections 1-3 (T1 recursion, initial data).
Initial integrand (sqrt part only; the AB part has no pole at v=i):
  Phi(v) = -(sqrt2/2) 8^d v^{2d-1} S_p(t(v)) / [(v^2-a^2)(v^2-b^2)(v^2+1)^{2d-2}]
  rho_{pq}(delta) = [s^delta] Common(s) * [S_p(t(i+s)) + [p>q] S_q(1-t(i+s))]
  Common(s) = -(sqrt2/2) 8^d (i+s)^{2d-1} (2i+s)^{2-2d} / G(i+s), G = v^4-6v^2+1.
  S_p(t(i+s)) = (1/(p!q!)) sum_j C_j sum_s Stir2(p,s) Stir2(q,j-s) M[(j,s)]
  M[(j,s)] = T^{2(j-s)} U^{2s} W^{d-j}, C_j = binom(1/2,j) j! 2^d,
  T = t(i+s), U = 1-T, W = q2(t(i+s)) = T^2+U^2.
"""
import os, sys
from fractions import Fraction as Fr
from math import comb, factorial

_HERE = os.path.dirname(os.path.abspath(__file__))
_WTS = os.path.normpath(os.path.join(_HERE, '..', 'WTS_weight_k_structure', 'theory'))
sys.path.insert(0, _WTS)
from fastred import QSI, ZERO, ONE, IQ, AQ, BQ

SQ2 = QSI(0,1,0,0)

# ---------- series arithmetic over QSI (series in s = v-i), truncated to s^D ----------
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
    res = [ONE] + [ZERO]*D
    if e == 0: return res
    b = list(base)
    while e:
        if e & 1: res = ser_mul(res, b, D)
        e >>= 1
        if e: b = ser_mul(b, b, D)
    return res

def ser_inv(A, D):
    a0i = A[0].inv()
    out = [a0i] + [ZERO]*D
    for n in range(1, D+1):
        acc = ZERO
        for k in range(1, n+1):
            if k < len(A) and not A[k].is_zero() and not out[n-k].is_zero():
                acc = acc + A[k]*out[n-k]
        out[n] = (ZERO - acc)*a0i if not acc.is_zero() else ZERO
    return out

# ---------- Stirling numbers of 2nd kind ----------
_STIR = None
def stirling_table(maxn):
    global _STIR
    if _STIR is not None and len(_STIR) > maxn:
        return _STIR
    S = [[0]*(maxn+1) for _ in range(maxn+1)]
    S[0][0] = 1
    for n in range(1, maxn+1):
        for k in range(1, n+1):
            S[n][k] = k*S[n-1][k] + S[n-1][k-1]
    _STIR = S
    return S

def binom_half(j):
    r = Fr(1)
    for t in range(j):
        r *= (Fr(1,2) - t)
    return r / factorial(j)

def keys_of_degree(d):
    return [(i, d-i) for i in range(d, (d-1)//2, -1)]

# ---------- level-d series data ----------
def level_series(d, D, S2):
    P = [QSI(-2,0,2,0), QSI(2,0,2,0), ONE] + [ZERO]*(D-2); P = P[:D+1]
    den = [QSI(0,0,4,0), QSI(4,0,0,0)] + [ZERO]*(D-1); den = den[:D+1]
    T = ser_mul(P, ser_inv(den, D), D)
    U = [ONE - T[0]] + [ZERO - c for c in T[1:]]
    W = [a+b for a,b in zip(ser_mul(T,T,D), ser_mul(U,U,D))]
    T2 = ser_mul(T, T, D); U2 = ser_mul(U, U, D)
    T2p = [[ONE]+[ZERO]*D]; U2p = [[ONE]+[ZERO]*D]; Wp = [[ONE]+[ZERO]*D]
    for a in range(1, d+1):
        T2p.append(ser_mul(T2p[-1], T2, D))
        U2p.append(ser_mul(U2p[-1], U2, D))
        Wp.append(ser_mul(Wp[-1], W, D))
    M = {}
    for j in range(1, d+1):
        Wdj = Wp[d-j]
        for s in range(0, j+1):
            M[(j,s)] = ser_mul(ser_mul(T2p[j-s], U2p[s], D), Wdj, D)
    C = [Fr(0)]*(d+1)
    for j in range(1, d+1):
        C[j] = binom_half(j) * factorial(j) * (2**d)
    ips = [IQ, ONE] + [ZERO]*(D-1); ips = ips[:D+1]
    v_pow = ser_pow_small(ips, 2*d-1, D)
    twoi = [QSI(0,0,2,0), ONE] + [ZERO]*(D-1); twoi = twoi[:D+1]
    twoi_neg = ser_pow_small(ser_inv(twoi, D), 2*d-2, D)
    v2 = ser_mul(ips, ips, D); v4 = ser_mul(v2, v2, D)
    G = [v4[t] - v2[t].scale(Fr(6)) + (ONE if t==0 else ZERO) for t in range(D+1)]
    Common = ser_mul(ser_mul(v_pow, twoi_neg, D), ser_inv(G, D), D)
    prefac = SQ2 * QSI(Fr(-(8**d), 2))
    Common = [c*prefac for c in Common]
    return dict(d=d, D=D, M=M, C=C, Common=Common)

def rho_key(LS, p, q, S2):
    D = LS['D']; d = LS['d']; M = LS['M']; C = LS['C']
    inv_pf = Fr(1, factorial(p)*factorial(q))
    Sp = [ZERO]*(D+1)
    for j in range(1, d+1):
        Cj = C[j]
        if Cj == 0: continue
        slo = max(0, j-q); shi = min(j, p)
        for s in range(slo, shi+1):
            Sps = S2[p][s]; Sqj = S2[q][j-s]
            if Sps == 0 or Sqj == 0: continue
            coeff = Cj * Sps * Sqj
            Ms = M[(j,s)]
            for t in range(D+1):
                if not Ms[t].is_zero():
                    Sp[t] = Sp[t] + Ms[t].scale(coeff)
    Sp = [c.scale(inv_pf) for c in Sp]
    out = ser_mul(LS['Common'], Sp, D)
    if p > q:
        Sq = [ZERO]*(D+1)
        for j in range(1, d+1):
            Cj = C[j]
            if Cj == 0: continue
            slo = max(0, j-p); shi = min(j, q)
            for s in range(slo, shi+1):
                Sqs = S2[q][s]; Spj = S2[p][j-s]
                if Sqs == 0 or Spj == 0: continue
                coeff = Cj * Sqs * Spj
                Ms = M[(j, j-s)]
                for t in range(D+1):
                    if not Ms[t].is_zero():
                        Sq[t] = Sq[t] + Ms[t].scale(coeff)
        Sq = [c.scale(inv_pf) for c in Sq]
        fold = ser_mul(LS['Common'], Sq, D)
        out = [a+b for a,b in zip(out, fold)]
    return out   # rho_{pq}(delta) = out[delta]

def rho_level(d, D, S2):
    LS = level_series(d, D, S2)
    return {key: rho_key(LS, key[0], key[1], S2) for key in keys_of_degree(d)}

# ---------- tau Taylor coefficients of DL1, DL2 at v=i ----------
def tau1(n):
    s = n+1
    return ((-1)**s) * ( (IQ-BQ)**(-(s+1)) + (IQ+AQ)**(-(s+1)) - IQ**(-(s+1)) )
def tau2(n):
    s = n+1
    return ((-1)**s) * ( (IQ-AQ)**(-(s+1)) + (IQ+BQ)**(-(s+1)) - IQ**(-(s+1)) )

# ---------- T1 recursion (towlib.Lambda_ij), depth-truncated ----------
def Lambda_ij(i, j, delta, U):
    """(Lambda_d U)(i,j) at depth delta; U(key, dp) = degree-(d+1) depth-dp data, d=i+j."""
    d = i+j
    tot = ZERO
    for dp in range(0, delta+1):
        c1 = tau1(delta-dp).scale(Fr(i+1, 2*d-dp-1))
        tot = tot + U((i+1,j), dp) * c1
        if i > j:
            c2 = tau2(delta-dp).scale(Fr(j+1, 2*d-dp-1))
            tot = tot + U((i,j+1), dp) * c2
    if j == i-1:
        m = 4*i-4-delta
        ft = ZERO
        for mu in range(m, 4*i-3):
            g = ZERO
            for dp in range(0, 4*i-4-mu+1):
                g = g - U((i,i), dp) * tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
            ft = ft + g.conj_i() * (IQ**(-mu)) * Fr(comb(mu-2, mu-m)) * (IQ**(m-2))
        tot = tot + ft.scale(-i)
    return tot

def build_tower_level(k, D, rho_k, S2=None):
    """tower at level k: dict key -> [u(0..D)] for all keys of degree 2..k.
       rho_k: dict degree-k key -> [rho(0..D)]."""
    tower = {key: list(vals) for key, vals in rho_k.items()}
    for n in range(k-1, 1, -1):
        for key in keys_of_degree(n):
            i, j = key
            def U(kk, dp, _t=tower):
                v = _t.get(kk)
                if v is None: return ZERO
                return v[dp] if dp < len(v) else ZERO
            tower[key] = [Lambda_ij(i, j, delta, U) for delta in range(D+1)]
    return tower
'''
with open('trunc_engine.py','w') as f:
    f.write(engine)
print("wrote trunc_engine.py", len(engine), "bytes")
