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
mp.mp.dps = 60

def hmp(w):
    hv = cc.h(w)
    return mp.mpc(mp.mpf(hv.real.mid().str(55,radius=False,more=True)), mp.mpf(hv.imag.mid().str(55,radius=False,more=True)))

Li1i = -mp.log(1-1j)
print('h((2,2))       :', mp.nstr(hmp((2,2)), 15))
print('Li1(i)^2/2     :', mp.nstr(Li1i**2/2, 15))
print('zeta(1,1;1,i)  :', mp.nstr(Li1i**2/2 - mp.polylog(2,-1)/2, 15), ' (formula claim)')

print()
print('h((4,4))       :', mp.nstr(hmp((4,4)), 15))
print('Li1(-1)^2/2 - Li2(1)/2:', mp.nstr(mp.log(2)**2/2 - mp.pi**2/12, 15))
print('Li1(-1)^2/2 (shuffle) :', mp.nstr(mp.log(2)**2/2, 15))
