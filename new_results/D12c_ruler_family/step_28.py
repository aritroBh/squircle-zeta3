import mpmath as mp
mp.mp.dps = 80

# H_q(eps) = ∫_0^1 u^{eps-1}(1+u)^{-1-eps} g_eps(u) du, g_eps(u) = (1+u^{q(1-eps)})^{1/q}-1
# K_m(q) = (-log2)^{m+1}/(m+1)! + H_q^{(m)}(0)/m!  ... wait: G = 2^{-eps} + eps*H(eps)
# G = sum (-log2)^k eps^k/k! + eps * sum H^{(j)}(0) eps^j/j!
# K_m = (-log2)^{m+1}/(m+1)! + H^{(m)}(0)/m!

def Hderiv(q, m, dps=80):
    """H^{(m)}(0) via Cauchy integral around eps=0 (radius 1/2)."""
    mp.mp.dps = dps
    q = mp.mpf(q)
    def H(e):
        e = mp.mpc(e)
        f = lambda u: u**(e-1)*(1+u)**(-1-e)*((1+u**(q*(1-e)))**(1/q)-1)
        return mp.quad(f, [0, mp.mpf('0.5'), 1], maxdegree=8)
    r = mp.mpf('0.4')
    return mp.diff(H, mp.mpf(0), n=m, method='quad', radius=r) if False else None

# mp.diff with method='quad' does Cauchy; but integrand of H is real quad. Simpler: expand symbolically.
# d/deps [u^{eps-1}(1+u)^{-1-eps} ((1+u^{q(1-eps)})^{1/q}-1)]
# = u^{-1}(1+u)^{-1} * d/deps[ exp(eps*(log u - log(1+u))) * (exp((1/q)log(1+u^{q(1-eps)}))-1) ]
# Do it with sympy series in eps.
import sympy as sp
eps_s, u_s, q_s = sp.symbols('eps u q', positive=True)

def Km_integrand(qval, m):
    """Return sympy expression in u for the coefficient: K_m(q) = (-log2)^{m+1}/(m+1)! + ∫_0^1 R_m(u)/(u(1+u)) du
    where R_m(u) = [eps^m] ( u^eps(1+u)^{-eps} ((1+u^{q(1-eps)})^{1/q}-1) ) ... times m!? Careful:
    eps*H(eps) = eps ∫ ... ; K_m from eps*H: coefficient of eps^m in eps*H is [eps^{m-1}]H; H = ∫ u^{-1}(1+u)^{-1} A(eps,u) du
    with A = u^eps (1+u)^{-eps} g_eps(u). So K_m = (-log2)^{m+1}/(m+1)! + ∫ [eps^{m-1}]A · du/(u(1+u))."""
    q = sp.Rational(qval)
    e = eps_s
    A = sp.exp(e*(sp.log(u_s)-sp.log(1+u_s))) * (sp.exp(sp.log(1+u_s**(q*(1-e)))/q) - 1)
    ser = sp.series(A, e, 0, m+1).removeO().expand()
    return sp.expand(ser).coeff(e, m-1)

# Verify K_1(2) = M_1
M1_val = mp.mpf('2.02095276965635835694690816772523888426901720579')/8
mp.mp.dps = 60
R1 = Km_integrand(2, 1)
print("R1(u) =", sp.simplify(R1))
f = sp.lambdify(u_s, R1, 'mpmath')
mp.mp.dps = 60
K1_2 = mp.log(2)**2/2 + mp.quad(lambda u: f(u)/(u*(1+u)), [0, mp.mpf('0.5'), 1], maxdegree=10)
print("K_1(2) =", mp.nstr(K1_2, 40))
print("M_1    =", mp.nstr(M1_val, 40))
print("|diff| =", mp.nstr(abs(K1_2-M1_val),3))
