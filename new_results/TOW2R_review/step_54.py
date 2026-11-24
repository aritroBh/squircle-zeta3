import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
os.chdir(os.path.join(workdir, 'tow2_copy'))
from fractions import Fraction as Fr
import genintegrand as GI
from fastred import QSI, ZERO, ONE, IQ, AQ, BQ
import towlib  # TOW2's copy: exact towers + rho from binteg pickles

def qsi_pmul(P, Q):
    out = [ZERO]*(len(P)+len(Q)-1)
    for a,ca in enumerate(P):
        if ca.is_zero(): continue
        for b,cb in enumerate(Q):
            if not cb.is_zero(): out[a+b] = out[a+b] + ca*cb
    while out and out[-1].is_zero(): out.pop()
    return out
def qsi_padd(P, Q):
    n = max(len(P), len(Q)); out = [ZERO]*n
    for k,c in enumerate(P): out[k]=out[k]+c
    for k,c in enumerate(Q): out[k]=out[k]+c
    return out
def qsi_pscale(P, q): return [c*q for c in P]
def qsi_peval_shift(P, v0):
    # return Q(y) = P(v0 + y) as QSI polynomial in y
    out = [ZERO]
    ypow = [ONE]
    for c in P:
        out = qsi_padd(out, qsi_pscale(ypow, c))
        ypow = qsi_pmul(ypow, [v0, ONE])
    return out

def Stilde(S):
    D = len(S)-1
    TV = [QSI(-1), QSI(2), QSI(1)]
    out = [ZERO]
    for s, c in enumerate(S):
        if c == 0: continue
        term = [ONE]
        for _ in range(s): term = qsi_pmul(term, TV)
        # (4v)^{D-s} = 4^{D-s} v^{D-s}
        term = qsi_pscale(term, QSI(Fr(4)**(D-s)))
        term = [ZERO]*(D-s) + term
        out = qsi_padd(out, term)
    return out, D

def rho_first_principles(p, q, delta_max=4):
    """u_{pq}^{(d)}(delta) for delta=0..delta_max from my own Phi construction."""
    d = p+q
    S = GI.S_poly(d, p)
    St, D = Stilde(S)
    assert D == 2*d
    # Phi = -(sqrt2/2) 8^d 4^{-D} v^{2d-1-D} Stilde (v^2+1)^{2-2d} / ((v^2-a^2)(v^2-b^2))
    c = QSI(0, Fr(-(8**d), 2*4**D))     # -(sqrt2/2) 8^d/4^D : sqrt2 * (...)
    num = qsi_pscale(St, c)
    # v^{2d-1-D} = v^{-1}: to denominator
    den = [QSI(0), ONE]   # v
    # (v^2-a^2)(v^2-b^2) = v^4 - 6v^2 + 1 (a^2+b^2=6, a^2 b^2=1)
    den = qsi_pmul(den, [ONE, ZERO, QSI(-6), ZERO, ONE])
    # (v^2+1)^{2d-2}
    v2p1 = [ONE, ZERO, ONE]
    for _ in range(2*d-2): den = qsi_pmul(den, v2p1)
    # shift v -> i + y
    Ny = qsi_peval_shift(num, IQ)
    Dy = qsi_peval_shift(den, IQ)
    # Dy has zero constant terms up to y^{2d-2}
    k0 = 0
    while k0 < len(Dy) and Dy[k0].is_zero(): k0 += 1
    assert k0 == 2*d-2, (k0, 2*d-2)
    D1 = Dy[k0:]
    # invert D1 mod y^(delta_max+1)
    n = delta_max+1
    inv = [ZERO]*n
    inv[0] = D1[0].inv()
    for k in range(1, n):
        acc = ZERO
        for j in range(1, k+1):
            if j < len(D1): acc = acc + D1[j]*inv[k-j]
        inv[k] = -acc * inv[0]  # 1/D1: inv[k] = -(sum_{j>=1} D1_j inv_{k-j})/D1_0
        inv[k] = inv[k] * ONE
    inv = [inv[k] / D1[0] if False else inv[k] for k in range(n)]
    # fix: recurrence should be inv[k] = -(sum D1_j inv_{k-j}) * (1/D1_0):
    inv = [ZERO]*n
    d0inv = D1[0].inv()
    inv[0] = d0inv
    for k in range(1, n):
        acc = ZERO
        for j in range(1, k+1):
            if j < len(D1): acc = acc + D1[j]*inv[k-j]
        inv[k] = -acc * d0inv
    # R = Ny * inv mod y^n
    R = [ZERO]*n
    for a in range(min(len(Ny), n)):
        for b in range(n-a):
            R[a+b] = R[a+b] + Ny[a]*inv[b]
    factor = Fr(2) if p > q else Fr(1)
    return [R[delta].scale(factor) for delta in range(n)]

bad = []; tested = 0
for dd in range(4, 11):
    for key in towlib.keys_of_degree(dd):
        p, q = key
        mine = rho_first_principles(p, q, 4)
        for delta in range(0, 5):
            if delta > 2*dd-4: continue
            ref = towlib.rho(key, delta)
            tested += 1
            if mine[delta] != ref:
                bad.append((key, delta, mine[delta], ref))
print(f'first-principles Phi -> rho(delta<=4) vs binteg data: {tested} cells (d=4..10, all keys), {len(bad)} bad')
for b in bad[:6]: print(b)
with open(os.path.join(workdir,'progress.log'),'a') as f:
    f.write(f'step_12: first-principles Phi_p extraction of rho(0..4) vs binteg: {tested} cells, {len(bad)} bad\n')
