"""
H1 (n = 3): curve of candidate Borel singularities tau(beta) = rho(beta) * zeta*(beta) for the ray family
l = (0, -beta L, -L), beta in [0,1]  (beta = 0: edge direction, beta = 1: vertex direction), rho = 1+beta,
zeta*(beta) = eta L at the first complex zero of e^{2eta m_1}+e^{2eta m_2}+e^{2eta m_3}.
With x = e^{3s}, s = 2 eta L/3:   x + x^{1-beta} + 1 = 0,   tau = (1+beta) (1/2) log x.
"""
import mpmath as mp
mp.mp.dps = 30
x = mp.mpc(-0.5, 1e-30)  # beta = 0 root, take the branch log x = log(1/2) + i pi
lx = mp.log(mp.mpf(1) / 2) + 1j * mp.pi
prev = lx
for i in range(0, 21):
    beta = mp.mpf(i) / 20
    # solve in y = log x : e^y + e^{(1-beta) y} + 1 = 0, continue from previous root
    y = mp.findroot(lambda y: mp.exp(y) + mp.exp((1 - beta) * y) + 1, prev)
    prev = y
    tau = (1 + beta) * y / 2
    print(f'beta={mp.nstr(beta, 3):>5}  log x = {mp.nstr(y, 10):>28}  tau = {mp.nstr(tau, 10):>28}  |tau| = {mp.nstr(abs(tau), 10)}')
