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

def S2(m):  # sum_{n1<n2<=m} 1/(n1 n2)
    return mp.nsum(lambda n2: mp.harmonic(n2-1)/n2, [2, m])

def term(n3):
    return (1j)**n3 * S2(n3-1)/n3
s = mp.nsum(term, [3, mp.inf], tol=mp.mpf('1e-35'))
print('series Li111(1,1,i):', mp.nstr(s, 20))
print('h((2,2,2))         :', mp.nstr(hmp((2,2,2)), 20))
print('diff:', mp.nstr(abs(s-hmp((2,2,2))), 4))
