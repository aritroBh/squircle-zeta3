import os, json
os.chdir("<project>/research/08_theory/AV6_squircle_lattice_zeta")
from mpmath import mp, mpf, nstr, zeta, log, power, dirichlet, pslq, gamma, pi
mp.dps = 60
# pull stored sweep values (files from earlier runs)
s2 = json.load(open("sweep2.json"))
print("stored sweep2 (p,s,Z,agree):")
for (p_, s_, Z, r) in s2:
    if s_ == 4:
        print(f"  p={p_}: Z_4 = {Z[:52]}  A-C={r[:8]}")
# PSLQ on special values
def beta(x): return dirichlet(x, [0,1,0,-1])
G_ = beta(2)
Z3 = mpf("6.8501698651610044331500297427010021769768944752926")
Zhalf = mpf("4.36271384859765549971325257061")
basis = [zeta(3), zeta(4), zeta(5), zeta(2)*G_, G_, log(2), log(2)**2, pi**3, pi**4, zeta(2)**2, gamma(mpf('1.25'))**4, gamma(mpf('1.1'))**10]
for name, Z in [("Z_3(4)", Z3), ("Z_{1/2}(4)", Zhalf)]:
    rel = pslq([Z]+basis)
    print(name, "PSLQ:", rel)
# also on Z_4(4) from sweep2
for (p_, s_, Z, r) in s2:
    if p_ == 4 and s_ == 4:
        Z4 = mpf(Z)
        print("Z_4(4) =", nstr(Z4, 52))
        print("Z_4(4) PSLQ:", pslq([Z4]+basis))
