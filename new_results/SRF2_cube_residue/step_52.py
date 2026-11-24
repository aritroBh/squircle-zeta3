import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from vlib import ts_nodes
from srf2_lib import expm1_s, tomp
ctx.dps = 80; mp.mp.dps = 80

def K2_small(eps, a0=mp.mpf('0.5'), Lam=mp.mpf('160')):
    """int_0^{a0} a^{2e+1} (1+a)^{-q-2} S0 da, termwise.  q=1+2e."""
    q = 1+2*eps
    s = -q-2
    Lmax = 220
    g1 = [mp.mpf(1)]
    for l in range(1, Lmax+1):
        g1.append(g1[-1]*(s-l+1)/l)
    bh = [mp.mpf(1)]
    for j in range(1, 60):
        bh.append(bh[-1]*(mp.mpf('0.5')-j+1)/j)
    two2e = 2-2*eps
    tot = mp.mpf(0)
    for j in range(0, 60):
        for l in range(0, Lmax+1):
            pw = l + two2e*j
            if pw > Lam: break
            c = g1[l]*bh[j]
            if c == 0: continue
            mu = 2*eps+1 + pw
            tot += c * a0**(mu+1)/(mu+1)
    return tot

def K2_big(eps, h=arb(1)/64, T=5.0):
    q = 1+2*eps
    a0 = arb('0.5')
    nodes = ts_nodes(h, T)
    s = arb(0)
    for (x, xm, w, par) in nodes:
        a = a0 + (1-a0)*x
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        s += (arb(2*eps+1)*la).exp()*(arb(-(q+2))*l1a).exp()*S0*w*(1-a0)
    return s

def K2_rel(eps_arb, eps_mp):
    return K2_small(eps_mp) + tomp(K2_big(eps_arb))

mp.mp.dps = 60
I2 = mp.quad(lambda a: (1+a)*a**3/(mp.sqrt(1+a**6)+1), [0, mp.mpf('0.3'), 1])
print("target K2(-2+) = -2 + I2' =", mp.nstr(-2+I2, 20), "  rho_{-2} =", mp.nstr(3*(-2+I2), 20))
for eta in ['0.1', '0.05', '0.02', '0.01', '0.005']:
    e_mp = -2 + mp.mpf(eta)
    K = K2_rel(arb(mp.nstr(e_mp,30)), e_mp)
    print(f"eta={eta}: K2 = {mp.nstr(K,20)}   3K2 = {mp.nstr(3*K,20)}")
