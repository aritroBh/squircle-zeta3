
# pwr1_rhombus.py -- PWR1 referee: exact rhombus engine from D11c's operator representation.
# P/(4r) = (1-tau) * sum_m g_{2m}(q) (-tau)^{2m}/(2m)! * B2m(tau),
# B2m = (d/dmu)^{2m} B(1-tau/2+mu, 1-tau/2-mu)|_0,  eta = -tau.
# log B(1-tau/2+mu,1-tau/2-mu) = sum_{k>=2} zeta_k [(-1)^k((mu-tau/2)^k+(-mu-tau/2)^k) - tau^k]/k - log(1-tau)
import sympy as sp

NMAX = 6
ZS = sp.symbols('Z2:'+str(NMAX+2))
Z = {k: ZS[k-2] for k in range(2, NMAX+2)}
q, mu, TAU = sp.symbols('q mu TAU')

def ser_dict(): return {}

def logBeta_exp():
    """exponent E(mu,tau) = sum_k zeta_k/k [(-1)^k((mu-tau/2)^k+(-mu-tau/2)^k) - tau^k], as poly in mu with ser-in-tau coeffs.
       Represent as dict {(m,k): coeff} meaning sum coeff * mu^m tau^k * (zeta factors)."""
    E = sp.Integer(0)
    for k in range(2, NMAX+1):
        Pk = sp.expand((-1)**k*((mu-TAU/2)**k + (-mu-TAU/2)**k) - TAU**k)
        E += Z[k]*Pk/k
    return E  # polynomial in mu, tau with Z coefficients, total tau-degree <= NMAX

def beta_series():
    """B(mu,tau) = exp(E) truncated to tau-degree NMAX, as polynomial in mu."""
    E = logBeta_exp()
    # truncate E to tau-degree NMAX
    def truncTau(e):
        e = sp.expand(e)
        return sum(c*TAU**i for i in range(NMAX+1) for c in [sp.expand(e).coeff(TAU, i)] if c != 0)
    E = truncTau(E)
    out = sp.Integer(1); term = sp.Integer(1)
    for m in range(1, NMAX+1):
        term = sp.expand(term*E)
        # drop mu-degree > NMAX and tau-degree > NMAX
        term = truncTau(term)
        term = sp.expand(term)
        term = sum(sp.expand(term).coeff(mu, j)*mu**j for j in range(0, NMAX+1) if sp.expand(term).coeff(mu, j) != 0)
        out = out + term/m  # wait: need factorial: handle below
        out = sp.expand(out)
    return out  # NOT USED (see beta_series2)

def beta_series2():
    E = logBeta_exp()
    def truncTau(e):
        e = sp.expand(e)
        return sum(sp.expand(e).coeff(TAU, i)*TAU**i for i in range(NMAX+1))
    E = truncTau(E)
    out = sp.Integer(1); term = sp.Integer(1)
    for m in range(1, NMAX+1):
        term = truncTau(sp.expand(term*E))
        term = sp.Rational(1, m)*term
        out = sp.expand(out + term)
        out = truncTau(out)
        out = sum(sp.expand(out).coeff(mu, j)*mu**j for j in range(0, NMAX+1))
    return sp.expand(out)

def g_even(qv):
    """g_e(z) = even part of sqrt(cosh z + q sinh z) = sqrt(cosh z) * sum_j binom(1/2,2j) (q tanh z)^{2j}; return coeffs g_{2m}(q): g_e = sum g_{2m} z^{2m}/(2m)!."""
    z = sp.symbols('z')
    cosh_s = sum(z**(2*i)/sp.factorial(2*i) for i in range(NMAX+2))
    sqrtcosh = sp.series(sp.sqrt(cosh_s), z, 0, NMAX+1).removeO().expand()
    tanh_s = sp.series(sp.tanh(z), z, 0, NMAX+1).removeO().expand()
    H = sum(sp.binomial(sp.Rational(1,2), 2*j) * (qv*tanh_s)**(2*j) for j in range(0, NMAX//2+1))
    ge = sp.expand(sqrtcosh * H)
    ge = sum(ge.coeff(z, i)*z**i for i in range(0, NMAX+1, 2))
    return {m: sp.expand(sp.factorial(2*m)*ge.coeff(z, 2*m)) for m in range(NMAX//2+1)}

def logP_rhombus():
    Bser = beta_series2()  # B(mu, tau), poly in mu, tau with Z coeffs
    g2 = g_even(q)
    tot = sp.Integer(0)
    for m in range(0, NMAX//2+1):
        B2m = sp.expand(sp.diff(Bser, mu, 2*m).subs(mu, 0))
        tot += g2[m] * (-TAU)**(2*m)/sp.factorial(2*m) * B2m
    # tot = P/(4r(1-tau))... then log(P/(4r)) = log(1-tau) + log(tot) -- but tot already includes -log(1-tau) in B? No:
    # B(mu,tau) as built = Gamma(1-tau/2+mu)Gamma(1-tau/2-mu)/Gamma(2-tau) * (1-tau)  [we dropped -log(1-tau)]
    # i.e. we used logB = zeta-part only, so B_built = B_true * (1-tau). Then P/(4r) = (1-tau)*int... = B_built-based sum.
    # Careful: P/(4r) = (1-tau) * sum g (-tau)^{2m} B^{(2m)}_true/(2m)!. B_true = B_built/(1-tau). So P/(4r) = sum g (-tau)^{2m} B_built^{(2m)}/(2m)!.
    # Hence log(P/(4r)) = log(tot). Expand:
    s = sp.series(sp.log(tot), TAU, 0, NMAX+1).removeO().expand()
    return s
