
from flint import acb, arb, ctx
import json
def I_cert(eta_str, U=160):
    ctx.prec = 260
    eta = acb(eta_str)
    def f(u, ana):
        return (1+(-u).exp())**(-2-eta) * 2 * (-eta*u/2).exp() * (eta*u).cosh().sqrt(analytic=ana) * (-u).exp()
    v = acb.integral(f, acb(0), acb(U), abs_tol=arb('1e-62'))
    em = abs(float(eta_str))
    # |K| <= 1 for eta >= -2; |H(eta u)| = 2 e^{-eta u/2}|sqrt cosh| <= 2 e^{|eta|u/2} e^{|eta|u/2} = 2 e^{|eta| u}
    tail = arb(2)*arb(-float(U)).exp()/arb(1-em)
    return v, tail
refs = {}
for eta in ['0.05','0.1','0.2','0.3','-0.05','-0.1','-0.2','-0.3']:
    v, tail = I_cert(eta)
    total_rad = arb(v.rad()) + tail
    refs[eta] = (str(v), str(total_rad))
    print(f"eta={eta}: I = {v}")
    print(f"          total radius (incl tail) = {total_rad}")
json.dump(refs, open("iref_certified.json","w"))
