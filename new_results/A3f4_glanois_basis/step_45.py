import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 1600
cc.N = 1000
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 50

def hmp(w):
    hv = cc.h(w)
    return mp.mpc(mp.mpf(hv.real.mid().str(45,radius=False,more=True)), mp.mpf(hv.imag.mid().str(45,radius=False,more=True)))

# series value of Li_{1,1}(1,i) = sum_{n2>=2} i^{n2} H_{n2-1}/n2, accelerated
s = mp.nsum(lambda n2: (1j)**n2 * mp.harmonic(n2-1)/n2, [2, mp.inf], tol=mp.mpf('1e-40'))
print('series Li_{1,1}(1,i):', mp.nstr(s, 20))
print('h((2,2))            :', mp.nstr(hmp((2,2)), 20))
print('diff:', mp.nstr(abs(s - hmp((2,2))), 4))
# and Li_{1,1}(i,i) for reference
s2 = mp.nsum(lambda n2: (1j)**n2 * mp.nsum(lambda n1: (1j)**n1/n1, [1, n2-1])/n2, [2, mp.inf], tol=mp.mpf('1e-30'))
Li1i = -mp.log(1-1j)
print('series Li_{1,1}(i,i):', mp.nstr(s2, 15), ' vs (Li1(i)^2-Li2(i^2))/2:', mp.nstr((Li1i**2-mp.polylog(2,-1))/2, 15))
