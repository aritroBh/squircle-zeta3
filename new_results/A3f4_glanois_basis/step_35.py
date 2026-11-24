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
mp.mp.dps = 50

# zeta(1,1,1; 1,1,i) = sum i^n3 (H^2 - H2)/(2 n3)
def direct(M=60000):
    s = mp.mpc(0)
    H = mp.mpf(0); H2 = mp.mpf(0)
    z = 1j
    for n in range(1, M+1):
        # H_{n-1}, H2_{n-1}
        s += z**n * (H*H - H2) / (2*n)
        H += 1/mp.mpf(n); H2 += 1/mp.mpf(n)**2
    return s  # tail O(H^2/n^2)

d = direct()
hv = cc.h((2,2,2))
hmp = mp.mpc(mp.mpf(hv.real.mid().str(45,radius=False,more=True)), mp.mpf(hv.imag.mid().str(45,radius=False,more=True)))
print('direct:', mp.nstr(d, 20))
print('h(2,2,2):', mp.nstr(hmp, 20))
print('diff:', mp.nstr(abs(d-hmp), 4))
