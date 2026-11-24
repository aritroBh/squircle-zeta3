from flint import arb, ctx
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt()
c_w = arb("1.0")
ecw = (rt2/pi).sqrt()*(-c_w/2).exp()
print("ecw =", ecw)
print("sqrt(2/pi) =", (rt2/pi).sqrt())
print("e^{-0.5} =", (-c_w/2).exp())
print("product check =", (rt2/pi).sqrt()*arb("0.6065306597126334236037995349911804534419181354871869556828921587"))
