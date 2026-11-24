from flint import acb, arb, ctx
ctx.prec = 300
out = acb(0)
for a, chi in [(1,1),(3,-1),(5,-1),(7,1)]:
    z = acb.zeta(acb(5), acb(a,8))
    print(a, chi, z.str(15))
    out += chi*z
print('sum', out.str(15))
print('div 8^5', (out/acb(8)**5).str(15))
