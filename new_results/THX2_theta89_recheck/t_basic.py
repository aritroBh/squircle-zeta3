from mzv import *
# basic sanity: DS/DUAL of small words, and numeric float check at weights<=5
import itertools, random, math
def allconv(w):
    return [x for x in itertools.product((0,1,-1), repeat=w) if conv(x)]
def val(w, N=20000):
    # numeric zeta-sum value of word via nested sums (float, low accuracy)
    import numpy as np
    S = word_to_sum(w)
    n = np.arange(1, N+1, dtype=float)
    # S[0] is outermost (largest n). innermost last.
    inner = None
    # compute from the last: A_d(n)=e^n/n^s ; cumulative sums strictly below
    cur = None
    for (s, e) in reversed(S):
        term = (float(e)**n if e == 1 else (-1.0)**n)/n**s
        if cur is not None: term = term*cur
        cs = np.cumsum(term); cur = np.concatenate(([0.0], cs[:-1]))  # sum over m<n
        last = cs
    return (-1)**len(S)*last[-1]*1.0
random.seed(1)
worst = 0
for W in (3, 4, 5):
    cw = allconv(W)
    for _ in range(40):
        u = random.choice(allconv(2) if W > 2 else cw); v = random.choice(allconv(W-2))
        r = DS(u, v)
        # use alternating partial-sum averaging crude: compare with N and 2N
        tot = sum(c*val(w) for w, c in r.items())
        worst = max(worst, abs(tot))
    for w in random.sample(cw, min(30, len(cw))):
        r = DUAL(w)
        tot = sum(c*val(x) for x, c in r.items()); worst = max(worst, abs(tot))
print("worst float residual of sampled DS/DUAL relations (N=20000, divergence-slow ones ~1e-3):", worst)
