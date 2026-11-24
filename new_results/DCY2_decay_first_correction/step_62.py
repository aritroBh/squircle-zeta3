from flint import arb, acb, ctx
ctx.dps = 40
z = acb("0.6")
w = (1-(-2*z).exp())/(2*z)
print("w =", w)
print("sqrt w =", w.sqrt())
z2 = acb(arb("0.6","0.7"))
w2 = (1-(-2*z2).exp())/(2*z2)
print("w2 =", w2)
print("sqrt w2 =", w2.sqrt())
# maybe the issue: acb(u) where u is arb interval
u = arb("0.6","0.7")
z3 = acb(u)
print("z3 =", z3)
print((1-(-2*z3).exp())/(2*z3))
