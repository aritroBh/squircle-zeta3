"""Evaluate level-8 word vectors numerically (h-word machinery from A3c check_conversion.py,
self-contained) and compare with quadrature."""
import sys, pickle, json, time
sys.dont_write_bytecode = True
from functools import cache
from fractions import Fraction as F
from collections import Counter
from flint import acb, arb, acb_series, acb_poly, ctx
import mpmath as mp

ctx.prec = 400
ctx.cap = 330
N = ctx.cap
mp.mp.dps = 100
roots = [acb(0, arb.pi()*j/4).exp() for j in range(8)]
roots[0], roots[2], roots[4], roots[6] = acb(1), acb(0, 1), acb(-1), acb(0, -1)

def af(c): return acb(str(c))

@cache
def sh(u, v):
    if not u: return ((v, 1),)
    if not v: return ((u, 1),)
    d = Counter()
    for w, m in sh(u[1:], v): d[(u[0],)+w] += m
    for w, m in sh(u, v[1:]): d[(v[0],)+w] += m
    return tuple(d.items())

@cache
def series(w, kind):
    if not w: return acb_series([1])
    assert w[-1] != -1
    tail = series(w[1:], kind)
    x = w[0]
    if x == -1:
        return acb_series([0]+[tail[j]/j for j in range(1, N)])
    pole = roots[x] if kind == 0 else (acb(1) if x == 8 else 1-roots[x])
    return (tail/acb_series([-pole, 1])).integral()

@cache
def half(w, kind):
    if not w: return acb(1)
    if w[-1] == -1:
        terms = dict(sh((-1,), w[:-1])); count = terms.pop(w)
        return (half(w[:-1], kind)*acb("0.5").log()
                - sum((a*half(t, kind) for t, a in terms.items()), acb(0)))/count
    return acb_poly(series(w, kind).coeffs())(acb("0.5"))

@cache
def g_at_one(w):
    out = acb(0)
    for k in range(len(w)+1):
        rev = tuple(-1 if x == 0 else 8 if x == -1 else x for x in reversed(w[:k]))
        out += (-1)**k*half(rev, 1)*half(w[k:], 0)
    return out

@cache
def h(w):
    return (-1)**sum(x != -1 for x in w)*g_at_one(tuple(-1 if x == -1 else (-x) % 8 for x in w))

def evaluate(d):
    return sum((af(c)*h(w) for w, c in d.items()), acb(0))

def main(which, quadkey):
    vec = pickle.load(open(f"wv_{which}.pkl", "rb"))
    vec = [{tuple(w): F(c) for w, c in d.items()} for d in vec]
    t0 = time.time()
    valr = evaluate(vec[0])
    print(f"{which} rational part evaluated ({time.time()-t0:.0f}s)", flush=True)
    vals = evaluate(vec[1])
    tot = valr + acb(2).sqrt()*vals
    q = json.load(open("quad_90.json"))
    expected = mp.mpf(q[quadkey])
    got = mp.mpf(tot.real.mid().str(90, radius=False, more=True))
    err = abs(got - expected)
    print(f"{which}: word-vector value = {mp.nstr(got, 60)}")
    print(f"   quad               = {mp.nstr(expected, 60)}")
    print(f"   |wv - quad| = {mp.nstr(err, 4)}   flint rad = {tot.real.rad().str(5)}")
    json.dump({"wv_value": mp.nstr(got, 90), "err_vs_quad": mp.nstr(err, 6)},
              open(f"wv_{which}_check.json", "w"))

main("w5", "W5")
