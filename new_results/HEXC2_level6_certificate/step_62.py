from flint import acb, acb_poly, arb, ctx
import inspect
print([m for m in dir(acb_poly) if not m.startswith('_')])
# test construction and ops
p = acb_poly([acb(1), acb(2), acb(3)])
print(p, p.length())
try:
    q = p.mullow(p, 5)
    print('mullow ok:', q)
except Exception as e:
    print('mullow err:', e)
try:
    print('p(2):', p(acb(2)))
except Exception as e:
    print('call err:', e)
print('coeff:', p[1], type(p[1]))
print('acb pi:', acb.pi())
print('exp:', (acb(0,1)*acb.pi()/3).exp())
