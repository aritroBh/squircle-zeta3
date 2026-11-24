import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from srf2_lib import expm1_s, tomp
ctx.dps = 80; mp.mp.dps = 80

def J_small(eps, a0=mp.mpf('0.5'), Lam=mp.mpf('160')):
    """int_0^{a0} a^{2e} [(1+a)^{-2-2e} S0 - 1] da via termwise integration (mpmath, exact)."""
    s = -2-2*eps
    # binomial coeffs binom(s, l) and binom(1/2, j)
    Lmax = 200
    g1 = [mp.mpf(1)]
    for l in range(1, Lmax+1):
        g1.append(g1[-1]*(s-l+1)/l)
    bh = [mp.mpf(1)]
    for j in range(1, 60):
        bh.append(bh[-1]*(mp.mpf('0.5')-j+1)/j)
    two2e = 2-2*eps
    tot = mp.mpf(0)
    # terms a^{l + two2e*j}, l>=0, j>=0, not both 0
    for j in range(0, 60):
        for l in range(0, Lmax+1):
            if l == 0 and j == 0: continue
            pw = l + two2e*j
            if pw > Lam: break
            c = g1[l]*bh[j]
            mu = 2*eps + pw
            tot += c * a0**(mu+1)/(mu+1)
    return tot

def J_big(eps, h=arb(1)/64, T=5.0, a0=arb('0.5')):
    """int_{a0}^1 a^{2e}[(1+a)^{-2-2e} S0 - 1] da via tanh-sinh on [a0,1]."""
    from vlib import ts_nodes
    nodes = ts_nodes(h, T)
    a0f = mp.mpf('0.5')
    s = arb(0)
    for (x, xm, w, par) in nodes:
        a = a0 + (1-a0)*x
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        f = S0*expm1_s(arb(-(2+2*eps))*l1a) + A2/(S0+1)
        s += (arb(2*eps)*la).exp()*f*w*(1-a0)
    return s

def J_reliable(eps_arb, eps_mp):
    J1 = J_small(eps_mp) + tomp(J_big(eps_arb))
    return 1/(2*eps_mp+1) + J1

print("check J(-1) = -1 + I = -0.84721308479397908661:")
print("J_reliable(-1) =", mp.nstr(J_reliable(arb(-1), mp.mpf(-1)), 25))
print()
for eta in ['0.2','0.1','0.05','0.02','0.01','0.005']:
    e_mp = -1 + mp.mpf(eta)
    e_arb = arb(e_mp.str(30))
    Jv = J_reliable(e_arb, e_mp)
    q = 1+2*e_mp
    print(f"eta={eta}: J={mp.nstr(Jv,20)}  -qJ (residue est) = {mp.nstr(-q*Jv, 20)}")
