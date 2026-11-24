import sympy as sp
p = sp.Symbol('p', integer=True, positive=True)
i = sp.I
t = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
def tc(side, n, conjg=False):
    if side == 1 and conjg: return (-1)**n * t[(2,n)]
    return t[(side,n)]
def rfun(dd, m, conjg=False):
    s = -2*dd if not conjg else 2*dd
    return 2*(dd-2*m) + s*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd, m, conjg=False):
    s = -1 if not conjg else 1
    return sp.Rational(1,2)*( (3*dd-3)/2 + s*i*rfun(dd,m,conjg)/2 - (dd-2)**2 )
# normalize by rho_pp(0) (diagonal, D=2p). Note conj(rho_pp(0)) = -rho_pp(0) (pure i*real).
# ratios of rho_{p+1,p}(0)/rho_pp(0):
rB0 = 2*i*(4*p-1)/(p+1)
def erho(delta, conjg=False):
    """rho_{p+1,p}(delta)/rho_pp(0); conjg applies i->-i to the RATIO (then overall sign from conj(rho_pp(0))=-rho_pp(0) handled separately)."""
    r0 = rB0 if not conjg else sp.conjugate(rB0)
    if delta == 0: return r0
    if delta == 1: return (-1 if not conjg else 1)*i*(2*p+1-2)*r0
    if delta == 2: return R2(2*p+1, p+1, conjg)*r0
def rho_pp(delta, conjg=False):
    """rho_pp(delta)/rho_pp(0)."""
    if delta == 0: return sp.Integer(1)
    if delta == 1: return (-1 if not conjg else 1)*i*(2*p-2)
    if delta == 2: return R2(2*p, p, conjg)
def ehat(delta, s, conjg=False):
    P = p+1
    if (delta,s)==(0,0): return erho(0,conjg)*tc(1,0,conjg)*P/(4*p-1)
    if (delta,s)==(1,0): return ( erho(1,conjg)*tc(1,0,conjg)*P/(4*p-2) + erho(0,conjg)*tc(1,1,conjg)*P/(4*p-1) )
    if (delta,s)==(2,0): return ( erho(2,conjg)*tc(1,0,conjg)*P/(4*p-3) + erho(1,conjg)*tc(1,1,conjg)*P/(4*p-2)
                                + erho(0,conjg)*tc(1,2,conjg)*P/(4*p-1) )
    if (delta,s)==(2,1): return (-erho(0,conjg)/(2*(2*p+1)-3))*tc(1,0,conjg)*P/(4*p-3)
def Delta(delta, conjg=False):
    """Delta_pp(delta)/rho_pp(0) = rho_pp(delta) - sum_s (-1)^s ehat(delta,s)."""
    v = rho_pp(delta, conjg) - sum( (-1)**s * ehat(delta, s, conjg) for s in range(0, delta//2+1) )
    return sp.expand(v)
def Delta_conj(delta):
    """conj(Delta_pp(delta))/rho_pp(0) = - conj_i(Delta/rho_pp(0)) since conj(rho_pp(0)) = -rho_pp(0)."""
    return -Delta(delta, conjg=True)
# feed into child (p,p-1) at depth delta: direct DL2 + fold
def feed(delta):
    # direct DL2: p * sum tau2_{delta-dp} Delta(dp)/(4p-3-dp)
    direct = p*sum( tc(2, delta-dp)*Delta(dp)/(4*p-3-dp) for dp in range(0, delta+1) )
    # fold: -p * sum_mu (-i)^mu binom(mu-2,mu-m) i^{m-2} conj(g_mu), m = 4p-4-delta
    m = 4*p-4-delta
    ft = 0
    for mu in range(m, 4*p-3):
        g = 0
        for dp in range(0, 4*p-4-mu+1):
            g = g - Delta_conj(dp)*tc(1, 4*p-4-mu-dp, True)*(-1)**(4*p-4-mu-dp)/(4*p-3-dp)
            # note: conj(g) uses conj(tau1_n) = (-1)^n tau2_n; tc(1,n,True) already = (-1)^n tau2_n
        g = sp.expand(g)
        ft = ft + g * (-i)**mu * sp.binomial(mu-2, mu-m) * i**(m-2)
    ft = -p*sp.expand(ft)
    return sp.expand(direct + ft)
for delta in (0,1,2):
    f = feed(delta)
    print(f'E3 delta={delta}: feed/rho_pp(0) =', sp.simplify(f))
