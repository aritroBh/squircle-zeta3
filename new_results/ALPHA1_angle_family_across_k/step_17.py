"""step_17.py -- T(r) piece by piece: T = -sum_{tau,rho} e_tau f_rho G_{rho,tau}(1),
G_{rho,tau}(1) = int_r^1 log((w-tau)/(r-tau))/(w-rho) dw  (regular at w=r).
Identify each G as dilog + log products."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 60

def pieces(r, verbose=False):
    e = [(0, mp.mpf(1)/2), (1, mp.mpf(1)/2), (-1, mp.mpf(-1)), (-r, mp.mpf(-1)/2)]
    f = [(r, mp.mpf(1)), (-1, mp.mpf(1)), (0, mp.mpf(-1))]
    out = {}
    tot = mp.mpf(0)
    for tau, et in e:
        for rho, fr in f:
            def g(w):
                return mp.log((w-tau)/(r-tau))/(w-rho)
            pts = [r, (r+1)/2, 1] if rho != r else [r, (3*r+1)/4, (r+1)/2, 1]
            G = mp.quad(g, pts)
            out[(tau, rho)] = G
            tot += -et*fr*G
    return out, tot

def T_alpha_c(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

# first verify the piece decomposition
for c in (mp.mpf(1)/2, mp.mpf(3)/5):
    r = (1-c)/(1+c)
    P, tot = pieces(r)
    T = T_alpha_c(c)
    print(f"c={mp.nstr(c,6)}: pieces sum={mp.nstr(tot,25)} T={mp.nstr(T,25)} diff={mp.nstr(abs(tot-T),3)}")

# identify each piece at several r
import itertools
rs = [mp.mpf(1)/3, mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(1)/9, mp.mpf(2)/5]
def cand_args(r, tau, rho):
    args = {}
    if tau != rho:
        args[f'Li2((1-rho)/(tau-rho))'] = (1-rho)/(tau-rho)
        args[f'Li2((r-rho)/(tau-rho))'] = (r-rho)/(tau-rho)
    return args
for (tau_r, rho_r) in [(0, r) for r in rs][:0]:
    pass
# For a fixed (tau, rho), print G at all r and PSLQ against [Li2 of the two args, pi^2, log products]
for tau_i, rho_i in [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2)]:
    Gs = []
    for r in rs:
        e = [(0, mp.mpf(1)/2), (1, mp.mpf(1)/2), (-1, mp.mpf(-1)), (-r, mp.mpf(-1)/2)]
        f = [(r, mp.mpf(1)), (-1, mp.mpf(1)), (0, mp.mpf(-1))]
        tau, et = e[tau_i]; rho, fr = f[rho_i]
        g = lambda w: mp.log((w-tau)/(r-tau))/(w-rho)
        pts = [r, (r+1)/2, 1] if rho != r else [r, (3*r+1)/4, (r+1)/2, 1]
        Gs.append(mp.quad(g, pts))
    # PSLQ at r = 1/3 with basis: Li2 args, pi^2, logs
    r = rs[0]
    e = [(0, mp.mpf(1)/2), (1, mp.mpf(1)/2), (-1, mp.mpf(-1)), (-r, mp.mpf(-1)/2)]
    f = [(r, mp.mpf(1)), (-1, mp.mpf(1)), (0, mp.mpf(-1))]
    tau, et = e[tau_i]; rho, fr = f[rho_i]
    B = [('1', mp.mpf(1))]
    for nm, a in cand_args(r, tau, rho).items():
        B.append((nm, mp.polylog(2, a)))
    # log products of endpoint diffs
    for x in {1-rho, r-rho, tau-rho, 1-tau, r-tau, abs(tau-rho), tau-rho}:
        B.append((f'log({mp.nstr(x,4)})', mp.log(x)))
    v = [Gs[0]] + [x for _, x in B]
    rel = mp.pslq(v, tol=mp.mpf(10)**(-50), maxcoeff=2000, maxsteps=3000)
    print(f"(tau={mp.nstr(tau,5)},rho={mp.nstr(rho,5)}): G={mp.nstr(Gs[0],20)}")
    if rel and rel[0] != 0:
        print("   rel:", [(str(rel[i]), B[i-1][0]) for i in range(1, len(rel)) if rel[i] != 0])
print("STEP17 DONE")
