import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 1200
cc.N = 800
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 40

def cmzv_word(xs, etas):
    p = len(xs); run = 0; signs=[0]*p
    for k in range(p-1,-1,-1):
        run = (run + etas[k]) % 8; signs[k]=run
    Tt=[]
    for k in range(p):
        Tt.append(signs[k]); Tt.extend([-1]*(xs[k]-1))
    return tuple(reversed(Tt))

i = 1j
def dzeta(xs, etas, M=30000):
    """direct nested sum for zeta(x; eps), eps = i^eta, i=exp(2 pi i/4)."""
    # eta as 8th-root exponents: zeta8 = exp(i pi/4); eps_k = zeta8^eta_k
    z = [mp.exp(mp.pi*mp.j/4*e) for e in etas]
    p = len(xs)
    if p == 1:
        return mp.nsum(lambda n: z[0]**n/n**xs[0], [1, mp.inf], tol=mp.mpf('1e-30'))
    if p == 2:
        # sum over n2, inner harmonic-like sums
        def outer(n2):
            inner = mp.nsum(lambda n1: z[0]**n1/n1**xs[0], [1, n2-1])
            return z[1]**n2 * inner / n2**xs[1]
        return mp.nsum(outer, [2, mp.inf], tol=mp.mpf('1e-25'))
    raise NotImplementedError

# depth-1 and depth-2 B4 elements
tests = [((3,),(2,)), ((1,2),(0,2)), ((2,1),(0,2))]
for xs, etas in tests:
    direct = dzeta(xs, etas)
    hv = cc.h(cmzv_word(xs, etas))
    hmp = mp.mpc(mp.mpf(hv.real.mid().str(35,radius=False,more=True)), mp.mpf(hv.imag.mid().str(35,radius=False,more=True)))
    print(xs, etas, 'diff:', mp.nstr(abs(direct-hmp),4), ' direct=', mp.nstr(direct,10), ' h=', mp.nstr(hmp,10))
