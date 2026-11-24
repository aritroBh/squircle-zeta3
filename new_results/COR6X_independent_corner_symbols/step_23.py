hpl = r'''
# hpl_num.py -- high-precision HPL evaluation H_w(x), letters in {0,1,-1}, |x|<1.
# Convention: H_{a,w}(x) = int_0^x f_a(t) H_w(t) dt (RV); the FIRST letter is integrated.
# Words ending in +-1 have a pure Taylor series (first-letter recursion).
# Words ending in 0 are reduced by exact shuffle (trailing zeros).
import mpmath as mp
from functools import lru_cache

mp.mp.dps = 80
GUARD = 40
KMAX = 6000

@lru_cache(maxsize=None)
def _series(word, prec):
    """Taylor coefficients c_k (mpf), k=1..KMAX, of H_word(x); word must end in +-1."""
    with mp.workprec(prec):
        word = tuple(word)
        if len(word) == 1:
            a = word[0]
            if a == 1:
                return tuple(mp.mpf(1) / k for k in range(1, KMAX + 1))
            if a == -1:
                return tuple(-mp.mpf((-1) ** k) / k for k in range(1, KMAX + 1))
            raise ValueError("bad word")
        a = word[0]
        c = _series(word[1:], prec)
        out = []
        if a == 0:
            for k in range(1, KMAX + 1):
                out.append(c[k - 1] / k)
        elif a == 1:
            s = mp.mpf(0)
            for k in range(1, KMAX + 1):
                if k >= 2:
                    s += c[k - 2]
                out.append(s / k)
        else:
            s = mp.mpf(0)
            for k in range(1, KMAX + 1):
                if k >= 2:
                    s = c[k - 2] - s
                out.append(s / k)
        return tuple(out)

def _eval_series(word, x, dps):
    cs = _series(tuple(word), mp.mp.dps + GUARD)
    with mp.workprec(dps + GUARD):
        x = mp.mpf(x)
        s = mp.mpf(0)
        xk = x
        for k, c in enumerate(cs, start=1):
            t = c * xk
            s += t
            xk *= x
            if k > 60 and abs(t) < mp.mpf(10) ** (-(dps + 15)):
                break
        return +s

def _shuffle(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + r for r in _shuffle(u[1:], v)] + [(v[0],) + r for r in _shuffle(u, v[1:])]

@lru_cache(maxsize=None)
def hpl(word, x, dps=None):
    word = tuple(word)
    if dps is None:
        dps = mp.mp.dps
    x = mp.mpf(x)
    if word[-1] != 0:
        return _eval_series(word, x, dps)
    # trailing zeros
    s = 0
    while s < len(word) and word[len(word) - 1 - s] == 0:
        s += 1
    if s == len(word):
        with mp.workprec(dps + GUARD):
            return +mp.log(x) ** s / mp.factorial(s)
    w = word[:-s]
    with mp.workprec(dps + GUARD):
        total = {}
        for ww in _shuffle(w, (0,) * s):
            total[ww] = total.get(ww, 0) + 1
        res = hpl(w, x, dps) * mp.log(x) ** s / mp.factorial(s)
        for ww, m in total.items():
            if ww == word and m == 1:
                continue
            res -= m * hpl(ww, x, dps)
        return +res
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/hpl_num.py","w").write(hpl)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_03.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("STDERR:", (r.stderr or "")[-2000:])
