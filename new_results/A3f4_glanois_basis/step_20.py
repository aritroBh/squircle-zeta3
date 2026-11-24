import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 2400
cc.N = 1200
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 40

xi8 = mp.exp(mp.pi*mp.j/4)

def z12_direct(eps1, eps2, M=40000):
    # sum_{n1<n2} eps1^n1 (eps2 xi8)^n2 / (n1 n2^2); eps=+-1
    z2 = eps2*xi8
    total = mp.mpc(0)
    # block over n2 mod 8 to accelerate: just plain sum to M, tail O(log M / M^2)
    H = mp.mpf(0)
    s = mp.mpc(0)
    for n2 in range(1, M+1):
        H += 1/mp.mpf(n2)   # H_{n2}
        if n2 >= 2:
            s += z2**n2 * (H - 1/mp.mpf(n2)) / n2**2
    # tail ~ z2^M (log M)/M^2
    return s

d = z12_direct(1, 1)
hw = cc.h((-1,1,1))
print('direct  :', mp.nstr(d, 20))
print('h(-1,1,1):', mp.nstr(cc.approx(hw), 20))
print('diff:', mp.nstr(abs(d-cc.approx(hw)),5))
