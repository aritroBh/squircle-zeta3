from flint import acb, ctx
import inspect
# find the right calling convention for Hurwitz zeta in python-flint
print([m for m in dir(acb) if 'zeta' in m.lower()])
ctx.prec = 300
try:
    v = acb.zeta(acb(5), acb(0.125))
    print('static ok', v.str(20))
except Exception as e:
    print('static fail', e)
