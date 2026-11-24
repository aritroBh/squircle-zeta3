from flint import arb, acb, ctx
ctx.dps = 30
a = arb(0.6, 0.7)
print("mid:", a.mid(), "rad:", a.rad())
b = arb(0.55, 0.05)
print("b:", b, "mid:", b.mid(), "rad:", b.rad())
# python-flint interval constructor:
c = arb(0.5, 0.6)   # ?
print("c:", c, c.mid(), c.rad())
# The proper way: arb(lo, hi) is mid, rad. For [lo,hi] use:
def iv(lo, hi):
    lo = arb(lo); hi = arb(hi)
    mid = (lo+hi)/2
    rad = (hi-lo)/2
    return arb(mid, rad)
d = iv("0.6","0.7")
print("iv:", d, d.mid(), d.rad())
print("acb:", acb(d))
