from flint import arb, acb, ctx
ctx.dps = 40
u = arb("0.6", "0.7")
def v_acb(z):
    return ((1-(-2*z).exp())/(2*z)).sqrt()
rho3 = v_acb(acb(u)) - 1 + u/2 - 5*u**2/24
print("rho3 =", rho3)
val = abs(rho3)/u**3
print("val =", val, type(val))
print("upper:", val.upper() if hasattr(val,'upper') else "n/a")
