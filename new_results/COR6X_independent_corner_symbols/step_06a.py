
# step_06a.py -- independent validation of hpl evaluator at x = 1/sqrt2 via Chebyshev integrals
import sys; sys.path.insert(0, '.')
import numpy as np
import mpmath as mp
import json, re
from hpl_num import hpl
mp.mp.dps = 50

N = 256
# Clenshaw-Curtis on [0, x0] via Chebyshev points
def cheb_grid(x0, n=N):
    k = np.arange(n)
    t = np.cos(np.pi * k / (n - 1))          # [-1,1]
    x = (t + 1) / 2 * x0                      # [0, x0]
    return x

def cheb_integrate(x, y):
    # integrate y(x) dx on grid x (descending Chebyshev from x0 to 0); use clenshaw-curtis weights
    n = len(x)
    # map to [-1,1]: x = (t+1)/2*x0; t from 1 to -1
    x0 = x[0]
    w = np.zeros(n)
    Nn = n - 1
    for j in range(n):
        # CC weights for t in [-1,1]
        s = 0.0
        for m in range(Nn // 2 + 1):
            cm = 1.0 if m == 0 else 2.0
            if 2 * m == Nn:
                cm *= 0.5
            s += cm * np.cos(2 * m * np.pi * j / Nn) / (1 - 4 * m * m)
        w[j] = (1.0 / Nn) * s * 2.0
    return 0.5 * x0 * float(np.dot(w, y))

x0 = float(1/np.sqrt(2))
grid = cheb_grid(x0)
from functools import lru_cache
@lru_cache(maxsize=None)
def G_grid(word):
    word = tuple(word)
    if len(word) == 1:
        a = word[0]
        if a == 0: return np.log(grid)
        return np.log((grid - a) / (0.0 - a))
    a = word[0]
    inner = G_grid(word[1:])
    y = inner / (grid - a)
    # cumulative integral int_0^x
    vals = np.array([cheb_integrate(grid[j::-1][::-1][: j + 1], y[j::-1][::-1][: j + 1]) if j >= 8 else
                     np.trapz(y[: j + 1], grid[: j + 1]) for j in range(len(grid))])
    return vals

# too slow per-point; do cumulative integration more cleverly:
def cum_integrate(grid, y):
    # returns F(x_j) = int_0^{x_j} y dx for each grid point, using CC on subgrids is O(N^2 log)
    n = len(grid)
    F = np.zeros(n)
    # use polynomial fit per interval of 8 points
    for j in range(1, n):
        # local integration from grid[j-1] to grid[j] with a few surrounding points (gauss via polyfit)
        i0 = max(0, j - 6); i1 = min(n, j + 2)
        xs = grid[i0:i1]; ys = y[i0:i1]
        # sort ascending
        order = np.argsort(xs); xs = xs[order]; ys = ys[order]
        p = np.polyfit(xs, ys, min(5, len(xs) - 1))
        P = np.polyint(p)
        F[j] = F[j - 1] + np.polyval(P, grid[j]) - np.polyval(P, grid[j - 1])
    return F

@lru_cache(maxsize=None)
def G_grid2(word):
    word = tuple(word)
    if len(word) == 1:
        a = word[0]
        if a == 0: return np.log(grid)
        return np.log((grid - a) / (0.0 - a))
    a = word[0]
    inner = G_grid2(word[1:])
    y = inner / (grid - a)
    return cum_integrate(grid, y)

data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]

c0 = 1/mp.sqrt(2)
bad = 0
for e, w, c in TERMS:
    hv = hpl(w, c0)
    gv = G_grid2(w)[-1]  # G at x0
    sgn = (-1)**sum(1 for t in w if t == 1)
    err = abs(float(hv) - sgn * gv)
    if err > 1e-10:
        bad += 1
        print("MISMATCH word", w, "hpl:", float(hv), "G*sgn:", sgn * gv, "err:", err)
print("words checked:", len(TERMS), " mismatches:", bad)
if bad == 0:
    O = mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)
    print("Omega(1/sqrt2) via validated evaluator =", mp.nstr(O, 50))
