from flint import arb, acb, ctx
ctx.dps = 50
x = arb('0.5')
# suspect: acb.log called on arb instance
try:
    y = acb(x).log()
    print("acb(x).log() =", y)
except Exception as e:
    print("exception:", e)
