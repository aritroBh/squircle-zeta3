hpl = r'''
# hpl_num.py -- high-precision HPL evaluation H_w(x), letters in {0,1,-1}, |x|<1.
# Exact-definition series recursion (mpf coefficients with guard digits),
# leading zeros removed by exact shuffle reduction.
import mpmath as mp
from functools import lru_cache

mp.mp.dps = 80
GUARD = 40

def _series_word(word):
    """Taylor coefficients c_k (mpf), k=1..K, of H_word(x) for word not starting with 0."""
    return _series_cached(tuple(word), mp.mp.dps + GUARD)

@lru_cache(maxsize=None)
def _series_cached(word, prec):
    with mp.workprec(prec):
        word = tuple(word)
        K = _K_needed(prec)
        if len(word) == 1:
            a = word[0]
            if a == 1:
                return tuple(mp.mpf(1) / k for k in range(1, K + 1))
            if a == -1:
                return tuple(-mp.mpf((-1) ** k) / k for k in range(1, K + 1))
            raise ValueError("single letter 0 is log x")
        a = word[-1]
        c = _series_cached(word[:-1], prec)
        if a == 0:
            return tuple(c[k - 1] / k for k in range(1, K + 1))
        out = []
        s = mp.mpf(0)
        for k in range(1, K + 1):
            # c'_k = (1/k) * sum_{j<k} c_j * (1 if a==1 else (-1)^{k-1-j})
            if a == 1:
                s = s + c[k - 2] if k >= 2 else s
                out.append(s / k)
            else:
                s = c[k - 2] - s if k >= 2 else mp.mpf(0)
                out.append(s / k)
        return tuple(out)

@lru_cache(maxsize=None)
def _K_needed(prec):
    # |c_k| = O(k^w); tail |x|^k: use q = 0.95 worst case handled by caller via actual x.
    return 4000

def _eval_series(word, x):
    cs = _series_word(word)
    with mp.workprec(mp.mp.dps + GUARD):
        x = mp.mpf(x)
        q = abs(x)
        s = mp.mpf(0)
        xk = x
        for k, c in enumerate(cs, start=1):
            s += c * xk
            xk *= x
            if k > 60 and abs(c * xk) < mp.mpf(10) ** (-(mp.mp.dps + 10)):
                break
        return +s

def _shuffle(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + r for r in _shuffle(u[1:], v)] + [(v[0],) + r for r in _shuffle(u, v[1:])]

@lru_cache(maxsize=None)
def hpl(word, x):
    """H_word(x) at mp.mp.dps precision (x mpf or float in (-1,1))."""
    word = tuple(word)
    x = mp.mpf(x)
    # split leading zeros
    r = 0
    while r < len(word) and word[r] == 0:
        r += 1
    if r == 0:
        return _eval_series(word, x)
    if r == len(word):
        with mp.workprec(mp.mp.dps + GUARD):
            return +mp.log(x) ** r / mp.factorial(r)
    u = word[r:]
    with mp.workprec(mp.mp.dps + GUARD):
        # shuffle(0^r, u) = (0^r u) + others with fewer leading zeros
        total = {}
        for w in _shuffle((0,) * r, u):
            total[w] = total.get(w, 0) + 1
        res = mp.log(x) ** r / mp.factorial(r) * hpl(u, x)
        for w, m in total.items():
            if w == word and m == 1:
                continue
            res -= m * hpl(w, x)
        return +res
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/hpl_num.py","w").write(hpl)
print("hpl_num.py written")
