import subprocess
code = r'''
"""rul1num.py -- independent numerical evaluators for level-2 MZVs.

E1 (rigorous, high precision): words over {0,-1} via Landen x->x/(x-1) to HPL series at 1/2
   (geometric rate 1/2; only the defining series are used).
E2: J_{a,k} = int_0^inf s^k log^a(1+e^{-s}) ds by 1D quadrature.
E3 (moderate precision): general convergent words via nested-series tails + Richardson/Shanks.
Named basis constants via mpmath.
"""
import mpmath as mp
from fractions import Fraction as Fr
import rul1

# ---------------- E1: words over {0,-1} via Landen ----------------

def _hpl_half_series(w, N):
    """H_w(1/2) for HPL word w over {0,1}: nested series sum (1/2)^{n_1}/prod n_i^{s_i}."""
    # parse blocks ending in 1
    S = []
    run = 0
    for b in w:
        if b == 0: run += 1
        else: S.append(run+1); run = 0
    assert run == 0, "HPL word must end in 1 for convergence at 1/2"
    mp.mp.dps = mp.mp.dps
    # tails: g[k][m] = sum_{m>n_k>...} prod_{j>=k} 1/n_j^{s_j}; g_last(m)=sum_{n=1}^{m}1/n^{s_d}
    d = len(S)
    g = [mp.mpf(0)]*(N+1)
    acc = mp.mpf(0)
    for n in range(1, N+1):
        acc += 1/mp.mpf(n)**S[-1]
        g[n] = acc
    for k in range(d-2, -1, -1):
        new = [mp.mpf(0)]*(N+1)
        acc = mp.mpf(0)
        sk = S[k]
        for n in range(1, N+1):
            acc += g[n-1]/mp.mpf(n)**sk if k > 0 else mp.mpf(0)
            new[n] = acc
        g = new
    # now g[n] = sum_{n>n_2>...} prod_{j>=2}... ; outer sum
    tot = mp.mpf(0)
    half = mp.mpf(1)/2
    xp = half
    s1 = S[0]
    for n in range(1, N+1):
        tot += xp*g[n-1]/mp.mpf(n)**s1 if d > 1 else xp/mp.mpf(n)**s1
        xp *= half
    return tot

def landen_transform(w):
    """H_w(x) at x = y/(y-1): dict w' -> eps (+1/-1), H_w(x) = sum eps H_{w'}(y)."""
    out = {(): 1}
    for letter in w:
        new = {}
        choices = [(1, -1)] if letter == 1 else [(0, 1), (1, 1)]
        for w0, e0 in out.items():
            for c, s in choices:
                nw = w0 + (c,)
                new[nw] = new.get(nw, 0) + e0*s
        out = new
    return {k: v for k, v in out.items() if v}

def g_word_over_0m1_num(w, N=None, dps=None):
    """G(w;1) for convergent word w over {0,-1}: = (-1)^{#(-1)} H_{w[-1->1]}(-1), Landen to 1/2."""
    if dps: mp.mp.dps = dps
    if N is None:
        N = int(mp.mp.dps*3.4) + 30
    hw = tuple(1 if a == -1 else a for a in w)
    tot = mp.mpf(0)
    for w2, e2 in landen_transform(hw).items():
        tot += e2*_hpl_half_series(w2, N)
    return (-1)**sum(1 for a in w if a == -1)*tot

# ---------------- E2: J-integrals ----------------

def J_num(a, k, dps=None):
    if dps: mp.mp.dps = dps
    f = lambda s: s**k*mp.log(1+mp.e**(-s))**a
    return mp.quad(f, [0, mp.inf])

# ---------------- E3: general word, moderate precision ----------------

def zeta_stuple_num(S, dps=40, K=13):
    """zeta(S) for convergent S via tail arrays + Richardson/Shanks. Moderate precision."""
    old = mp.mp.dps
    mp.mp.dps = dps + 15
    S = tuple(S)
    d = len(S)
    if d == 0:
        mp.mp.dps = old; return mp.mpf(1)
    N = 2**K
    eps = [e for (s, e) in S]
    ss = [s for (s, e) in S]
    # tails from the bottom: g[m] = sum_{m>n_d... }
    g = [mp.mpf(0)]*(N+1)
    acc = mp.mpf(0)
    ed, sd = eps[-1], ss[-1]
    for n in range(1, N+1):
        acc += (mp.mpf(ed)**n)/mp.mpf(n)**sd
        g[n] = acc
    for k in range(d-2, 0, -1):
        ek, sk = eps[k], ss[k]
        new = [mp.mpf(0)]*(N+1)
        acc = mp.mpf(0)
        for n in range(1, N+1):
            acc += (mp.mpf(ek)**n)*g[n-1]/mp.mpf(n)**sk
            new[n] = acc
        g = new
    # outer partial sums at 2^j
    e1, s1 = eps[0], ss[0]
    partials = []
    acc = mp.mpf(0)
    n = 1
    targets = [2**j for j in range(3, K+1)]
    ti = 0
    for n in range(1, N+1):
        acc += (mp.mpf(e1)**n)*g[n-1]/mp.mpf(n)**s1
        if ti < len(targets) and n == targets[ti]:
            partials.append(acc); ti += 1
    if e1 == 1:
        val, err = mp.richardson(partials)
    else:
        # Euler-average consecutive partial sums to damp oscillation, then richardson
        avg = [(partials[i]+partials[i+1])/2 for i in range(len(partials)-1)]
        val, err = mp.richardson(avg)
    mp.mp.dps = old
    return +val

def word_num(w, dps=40):
    """G(w;1) = (-1)^d zeta(stuple(w)) for convergent words (moderate precision)."""
    S = rul1.word_to_stuple(w)
    return (-1)**len(S)*zeta_stuple_num(S, dps=dps)

def combo_num(combo, dps=40):
    return sum(mp.mpf(c.numerator)/mp.mpf(c.denominator)*word_num(w, dps) for w, c in combo.items())

# ---------------- named basis constants (mpmath) ----------------

def named_values(dps):
    mp.mp.dps = dps
    l2 = mp.log(2); pi2 = mp.pi**2; pi4 = mp.pi**4
    z3 = mp.zeta(3); z5 = mp.zeta(5)
    Li4 = mp.polylog(4, mp.mpf(1)/2); Li5 = mp.polylog(5, mp.mpf(1)/2)
    return {
        '1': mp.mpf(1),
        'log2': l2,
        'pi2': pi2, 'log2^2': l2**2,
        'z3': z3, 'pi2*log2': pi2*l2, 'log2^3': l2**3,
        'pi4': pi4, 'z3*log2': z3*l2, 'pi2*log2^2': pi2*l2**2, 'log2^4': l2**4, 'Li4h': Li4,
        'z5': z5, 'z3*pi2': z3*pi2, 'pi4*log2': pi4*l2, 'z3*log2^2': z3*l2**2,
        'pi2*log2^3': pi2*l2**3, 'log2^5': l2**5, 'Li5h': Li5, 'Li4h*log2': Li4*l2,
    }
'''
import pathlib
pathlib.Path('rul1num.py').write_text(code)
print("written")
