import sympy as sp
w = sp.symbols('w')
z = sp.symbols('z', positive=True)
# log phi(w) = 3 log w + 2 logGamma(w) - logGamma(2w+1) + log(psi(2w+1)-psi(w))
# Expand each piece as series in w at 0
N = 14
lg2 = 2*sp.loggamma(w).series(w, 0, N+1).removeO()
lg3 = sp.loggamma(2*w+1).series(w, 0, N+1).removeO()
dpsi = (sp.polygamma(0, 2*w+1) - sp.polygamma(0, w))
# psi(w) has -1/w pole; handle: psi(2w+1)-psi(w) = -1/w + (regular)
# get series via manual: psi(w) = -1/w - gamma + sum zeta(k+1)(-w)^k ...
gam = sp.EulerGamma
psi_w = -1/w - gam + sum((-1)**(k+1)*sp.zeta(k+1)*w**k for k in range(1, N+2))
psi_2w1 = -gam + sum((-1)**k*sp.zeta(k+1)*(2*w)**k for k in range(1, N+2))  # psi(1+x)=-gam+sum(-1)^k zeta(k+1) x^k? check sign
# psi(1+x) = -gamma + sum_{k>=1} (-1)^{k+1} zeta(k+1) x^k
psi_2w1 = -gam + sum((-1)**(k+1)*sp.zeta(k+1)*(2*w)**k for k in range(1, N+2))
d = sp.series(psi_2w1 - psi_w, w, 0, N).removeO()
logphi = 3*sp.log(w) + lg2 - lg3 + sp.log(d)
logphi_ser = sp.series(logphi, w, 0, N).removeO()
phi_ser = sp.series(sp.exp(logphi_ser), w, 0, N).removeO().expand()
print(phi_ser)
# extract coefficients rationally
coeffs = {}
for k in range(0, N):
    c = sp.expand(phi_ser).coeff(w, k)
    coeffs[k] = sp.nsimplify(c) if c != 0 else 0
    print(k, sp.simplify(c))
