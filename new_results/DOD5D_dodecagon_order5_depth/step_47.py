import flint, time
flint.ctx.prec = 450
x = flint.arb("0.373737")
t0=time.time()
s = flint.arb(0)
for i in range(20000):
    s += x.log()*x
print("arb mul+log add: %.1f us/op" % ((time.time()-t0)/20000*1e6), s)
t0=time.time()
a = flint.arb("1.234"); b = flint.arb("-0.987")
s = flint.arb(0)
for i in range(100000):
    s += a*b
print("arb mul+add: %.2f us/op" % ((time.time()-t0)/100000*1e6))
