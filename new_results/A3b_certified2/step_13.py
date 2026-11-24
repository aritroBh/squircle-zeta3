
from flint import arb, acb, ctx
ctx.prec = 120

# 1) Does the draft's P() crash on the j=0 term? (mul() with no args)
def mul(*xs):
    out = xs[0]
    for x in xs[1:]:
        out = out * x
    return out
try:
    mul()
    print("mul() with no args: OK??")
except Exception as e:
    print("mul() with no args ->", type(e).__name__, e)

# 2) arb methods we need
x = arb(3) / arb(7)
print("rad:", x.rad(), "| upper:", x.upper(), "| abs_upper:", x.abs_upper())
z = acb(arb(1,2), arb(-1,1))
print("acb real rad:", z.real.rad(), "imag rad:", z.imag.rad())

# 3) analytic log behavior on a cut-crossing ball
zball = acb(arb(-1, 0.5))   # straddles negative real axis
print("log no flag :", zball.log())
print("log analytic:", zball.log(analytic=True))
print("is_finite analytic:", zball.log(analytic=True).is_finite())

# 4) acb.integral doc
help(acb.integral)
