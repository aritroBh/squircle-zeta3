import sympy as sp
p = sp.Symbol('p', integer=True, positive=True)
i = sp.I
t = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
def tc1c(n):  # conj(tau1_n) = (-1)^n tau2_n
    return (-1)**n * t[(2,n)]
def rfun(dd, m, conjg=False):
    s = -2*dd if not conjg else 2*dd
    return 2*(dd-2*m) + s*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd, m, conjg=False):
    s = -1 if not conjg else 1
    return sp.Rational(1,2)*( (3*dd-3)/2 + s*i*rfun(dd,m,conjg)/2 - (dd-2)**2 )
rB0 = 2*i*(4*p-1)/(p+1)
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
    if (delta,s)==(0,0): return erho(0,conjg)*tc1c(0)*0 + erho(0,conjg)*t[(1,0)]*(0 if conjg else 1)*P/(4*p-1) if False else erho(0,conjg)*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-1)
    if (delta,s)==(1,0): return ( erho(1,conjg)*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-2) + erho(0,conjg)*( (-t[(2,1)]) if conjg else t[(1,1)])*P/(4*p-1) )
    if (delta,s)==(2,0): return ( erho(2,conjg)*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-3)
                                + erho(1,conjg)*( (-t[(2,1)]) if conjg else t[(1,1)])*P/(4*p-2)
                                + erho(0,conjg)*( t[(2,2)] if conjg else t[(1,2)])*P/(4*p-1) )
    if (delta,s)==(2,1): return (-erho(0,conjg)/(4*p-1))*(t[(2,0)] if conjg else t[(1,0)])*P/(4*p-3)
def Delta(delta, conjg=False):
    return sp.expand( rho_pp(delta, conjg) - sum( (-1)**s * ehat(delta, s, conjg) for s in range(0, delta//2+1) ) )
def Dc(delta):  # conj(Delta)/rho_pp(0) = -conj_i(Delta/rho)
    return -Delta(delta, conjg=True)
def gfun(mu, delta):
    """conj(g_mu) for fold at depth delta (m = 4p-4-delta), per formulas."""
    m = 4*p-4-delta
    # g_mu = - sum_dp Delta(dp) tau1_{4p-4-mu-dp}/(4p-3-dp) ; conjugated:
    terms = {4*p-4: [(0,0)], 4*p-5: [(0,1),(1,0)], 4*p-6: [(0,2),(1,1),(2,0)]}
    tot = 0
    for (dp, n) in terms[mu]:
        tot = tot - Dc(dp)*tc1c(n)/(4*p-3-dp)
    return tot
def feed(delta):
    direct = p*sum( t[(2, delta-dp)]*Delta(dp)/(4*p-3-dp) for dp in range(0, delta+1) )
    m = 4*p-4-delta
    if delta == 0:
        ft = -p*( gfun(4*p-4, delta)*(-1) )
    elif delta == 1:
        ft = -p*( gfun(4*p-5, delta)*(-1) + gfun(4*p-4, delta)*i*(4*p-6) )
    elif delta == 2:
        ft = -p*( gfun(4*p-6, delta)*(-1) + gfun(4*p-5, delta)*i*(4*p-7) + gfun(4*p-4, delta)*sp.binomial(4*p-6,2) )
    return sp.expand(direct + ft)
for delta in (0,1,2):
    f = feed(delta)
    f = f.replace(sp.binomial(4*p-6,2), (4*p-6)*(4*p-7)/2)
    print(f'E3 delta={delta}: feed =', sp.simplify(f))
