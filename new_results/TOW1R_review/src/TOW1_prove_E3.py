from proof_common import *
p = sp.Symbol('p', integer=True, positive=True)
# E3 (diagonal blindness), depths 0,1,2; normalized by rho_pp(0). conj(rho_pp(0)) = -rho_pp(0).
rB0 = 2*i*(4*p-1)/(p+1)                      # rho_{p+1,p}(0)/rho_pp(0)
def erho(delta, conjg=False):
    r0 = rB0 if not conjg else sp.conjugate(rB0)
    if delta == 0: return r0
    if delta == 1: return (-1 if not conjg else 1)*i*(2*p-1)*r0
    if delta == 2: return R2(2*p+1, p+1, conjg)*r0
def rho_pp(delta, conjg=False):
    if delta == 0: return sp.Integer(1)
    if delta == 1: return (-1 if not conjg else 1)*i*(2*p-2)
    if delta == 2: return R2(2*p, p, conjg)
def ehat(delta, s, conjg=False):
    P = p+1
    if (delta,s)==(0,0): return erho(0,conjg)*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-1)
    if (delta,s)==(1,0): return ( erho(1,conjg)*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-2)
                                + erho(0,conjg)*((-t[(2,1)]) if conjg else t[(1,1)])*P/(4*p-1) )
    if (delta,s)==(2,0): return ( erho(2,conjg)*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-3)
                                + erho(1,conjg)*((-t[(2,1)]) if conjg else t[(1,1)])*P/(4*p-2)
                                + erho(0,conjg)*(t[(2,2)] if conjg else t[(1,2)])*P/(4*p-1) )
    if (delta,s)==(2,1): return (-erho(0,conjg)/(4*p-1))*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-3)
def Delta(delta, conjg=False):
    return sp.expand( rho_pp(delta, conjg) - sum( (-1)**s * ehat(delta, s, conjg) for s in range(0, delta//2+1) ) )
def Dc(delta):
    return -Delta(delta, conjg=True)
def gfun(mu):
    terms = {4*p-4: [(0,0)], 4*p-5: [(0,1),(1,0)], 4*p-6: [(0,2),(1,1),(2,0)]}
    return sum( -Dc(dp)*tc1c(n)/(4*p-3-dp) for (dp, n) in terms[mu] )
def feed(delta):
    direct = p*sum( t[(2, delta-dp)]*Delta(dp)/(4*p-3-dp) for dp in range(0, delta+1) )
    if delta == 0:
        ft = -p*( gfun(4*p-4)*(-1) )
    elif delta == 1:
        ft = -p*( gfun(4*p-5)*(-1) + gfun(4*p-4)*i*(4*p-6) )
    elif delta == 2:
        ft = -p*( gfun(4*p-6)*(-1) + gfun(4*p-5)*i*(4*p-7) + gfun(4*p-4)*sp.binomial(4*p-6,2) )
    return sp.expand(direct + ft)
for delta in (0,1,2):
    f = feed(delta).replace(sp.binomial(4*p-6,2), (4*p-6)*(4*p-7)/2)
    print(f'E3 delta={delta}: feed remainder =', sp.simplify(f))
