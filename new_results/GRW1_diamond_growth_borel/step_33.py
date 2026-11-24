
from flint import acb, arb, ctx, showgood
import mpmath as mp
mp.mp.dps = 60

def I_cert(eta_str, U=64):
    ctx.prec = 200
    eta = acb(eta_str)
    def f(u, ana):
        return (1+(-u).exp())**(-2-eta) * 2 * (-eta*u/2).exp() * (eta*u).cosh().sqrt(analytic=ana) * (-u).exp()
    v = acb.integral(f, acb(0), acb(U), abs_tol=arb('1e-45'))
    # rigorous tail bound on [U, inf): |integrand| <= 2*max(1,2^{|eta+2|}) * e^{-u}(1+|eta|)... use |H|<=2 e^{|eta| u}, |K| <= max(1,2^{|2+eta|})
    em = abs(float(eta_str))
    tail = arb(2)*arb(2)**arb(abs(2+em))*arb(-float(U)).exp()/arb(1-em)
    return v, tail

refs = {}
for eta in ['0.1','0.2','0.3','-0.1','-0.2','-0.3']:
    v, tail = I_cert(eta)
    refs[eta] = v
    print(f"eta={eta}: I = {v}")
    print(f"          rigorous tail <= {tail}")
import json
json.dump({k: str(v) for k,v in refs.items()}, open("iref_certified.json","w"))
