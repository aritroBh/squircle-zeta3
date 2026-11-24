"""Self-contained level-8 H-word evaluator (arb), ported from A3f4/Q/check_conversion.py.
IMPORTANT (A3f4 pitfall): set ctx.prec AND ctx.cap, then compute roots."""
from functools import cache
from collections import Counter
from flint import acb, arb, acb_series, acb_poly, ctx

ctx.threads = 1
ctx.prec = 500      # ~150 digits
ctx.cap = 900
N = ctx.cap
roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
roots[0],roots[2],roots[4],roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)

@cache
def sh(u, v):
    if not u: return {v:1}
    if not v: return {u:1}
    out = Counter()
    for w,a in sh(u[1:],v).items(): out[(u[0],)+w] += a
    for w,a in sh(u,v[1:]).items(): out[(v[0],)+w] += a
    return dict(out)

@cache
def series(w,kind):
    if not w: return acb_series([1])
    assert w[-1] != -1
    tail = series(w[1:],kind)
    x = w[0]
    if x == -1:
        return acb_series([0]+[tail[j]/j for j in range(1,N)])
    pole = roots[x] if kind == 0 else (acb(1) if x == 8 else 1-roots[x])
    return (tail/acb_series([-pole,1])).integral()

@cache
def half(w,kind):
    if not w: return acb(1)
    if w[-1] == -1:
        terms = dict(sh((-1,),w[:-1])); count = terms.pop(w)
        return (half(w[:-1],kind)*acb('0.5').log()
                -sum((a*half(t,kind) for t,a in terms.items()),acb(0)))/count
    return acb_poly(series(w,kind).coeffs())(acb('0.5'))

@cache
def g_at_one(w):
    out = acb(0)
    for k in range(len(w)+1):
        rev = tuple(-1 if x == 0 else 8 if x == -1 else x for x in reversed(w[:k]))
        out += (-1)**k*half(rev,1)*half(w[k:],0)
    return out

@cache
def h(w):
    w = tuple(w)
    return (-1)**sum(x != -1 for x in w)*g_at_one(tuple(-1 if x == -1 else (-x)%8 for x in w))

def evaluate_pairs(pairs):
    """pairs: iterable of (word, python int/fraction coeff) -> acb sum."""
    out = acb(0)
    for w,c in pairs:
        out += acb(str(c))*h(tuple(w))
    return out
