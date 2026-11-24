from mpmath import mp, mpf, zeta, log, pi, nprint, polylog, clsin, pslq, catalan
mp.dps = 210
z3 = zeta(3); L2 = log(2); Cl = clsin(2, pi/3); G = catalan
Li3_half = polylog(3, mpf('0.5')); Li2_half = polylog(2, mpf('0.5'))
odd = mpf(open('odd2.txt').read())
S1 = mpf('0.0502870873772545765534176332325716727641098374883945464791628')
# recompute S1' at 210 for precision
from mpmath import binomial
K = 800
S1p = mpf(0); Hk = mpf(0)
for k in range(0, K):
    if k>0: Hk += 1/mpf(k)
    H2k = Hk + sum(1/mpf(k+j) for j in range(1, k+1)) if k>0 else mpf(0)
    S1p += binomial(2*k,k)/mpf(16)**k/(2*k+1)*H2k
basis = [1, pi, pi**2, pi**3, z3, pi*z3, L2, pi*L2, pi**2*L2, L2**2, pi*L2**2, L2**3,
         Cl, pi*Cl, G, pi*G, Li3_half, Li2_half, pi*Li2_half, z3*L2, Cl*L2, log(3), pi*log(3),
         polylog(3, mpf('0.25'))]
names = ['1','pi','pi2','pi3','z3','pi*z3','log2','pi*log2','pi2*log2','log2^2','pi*log2^2',
         'log2^3','Cl2(pi/3)','pi*Cl2','G','pi*G','Li3(1/2)','Li2(1/2)','pi*Li2(1/2)',
         'z3*log2','Cl2*log2','log3','pi*log3','Li3(1/4)']

def hunt(target, label):
    v = [target] + basis
    rel = pslq(v, tol=mpf('1e-185'), maxcoeff=10**24, maxsteps=10000)
    print(f"--- PSLQ {label} ---")
    if rel is None:
        print("  no relation found")
    else:
        err = abs(sum(mpf(c)*x for c,x in zip(rel, v)))
        print("  residual:", nprint(err, 4))
        out = []
        if rel[0]: out.append(f"{rel[0]}*TARGET")
        for c,n in zip(rel[1:], names):
            if c: out.append(f"{c}*{n}")
        print("  " + " + ".join(out))
    return rel

hunt(odd, "odd part")
hunt(S1, "S1")
hunt(S1p, "S1'")
hunt(S1-S1p, "S1-S1'")
with open('progress.log','a') as fh:
    fh.write("step_24: PSLQ hunts on odd part, S1, S1', S1-S1' (24-dim basis, 200d)\n")
