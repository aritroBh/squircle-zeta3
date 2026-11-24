import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 4200
cc.N = 1400
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 50

def hmp(w, digs=45):
    hv = cc.h(w)
    return mp.mpc(mp.mpf(hv.real.mid().str(digs,radius=False,more=True)), mp.mpf(hv.imag.mid().str(digs,radius=False,more=True)))

# accelerated series for Li_{1,1,1}(1,1,i) = sum i^n3 (H^2-H2)/(2 n3)
def term(n3):
    H = mp.harmonic(n3-1); H2 = mp.harmonic(n3-1, 2)
    return (1j)**n3 * (H*H - H2)/(2*n3)
s = mp.nsum(term, [3, mp.inf], tol=mp.mpf('1e-40'))
print('series Li111(1,1,i):', mp.nstr(s, 20))
print('h((2,2,2))         :', mp.nstr(hmp((2,2,2)), 20))
print('diff:', mp.nstr(abs(s-hmp((2,2,2))), 4))
