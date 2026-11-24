
from flint import ctx, arb, acb
ctx.dps = 40
OMEGA = acb(1,1)/arb(2).sqrt()
N = arb(404, 4)
PI = arb.pi()
y = (PI/(2*N)).sqrt()
s0 = OMEGA/y

# H via the integral representation with interval s0
def H_int(s):
    def f(u, analytic):
        return (1-(-2*u).exp()).sqrt() * (-s*u).exp()
    return acb.integral(f, 0, 1) + acb.integral(f, 1, 10)
v = H_int(s0)
print("H(s0) integral, interval:", v.str(8))
print("radius:", v.real.rad(), v.imag.rad())

# H via gamma with interval s0
Hg = PI.sqrt()/4 * (s0/2).gamma() * (s0/2+arb(3)/2).rgamma()
print("H(s0) gamma, interval:", Hg.str(8))
print("radius:", Hg.real.rad(), Hg.imag.rad())

# digamma interval
d1 = (s0/2).digamma()
print("digamma interval:", d1.str(8))
